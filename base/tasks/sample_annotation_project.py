from __future__ import absolute_import, unicode_literals
import socket
from datetime import datetime
from celery import Task, shared_task, group, chain
from django.core.cache import cache
from django.conf import settings
import numpy as np
from base.models import SampleAnnotationProject, Molecule
from base.modules.cache_management import model_from_cache
from base.tasks.project import log_project, log_error
from fragmentation.utils import AdductManager

log_queue = settings.CELERY_LOG_QUEUE


class MetWorkTask(Task):
    def on_failure(self, exc, task_id, args, kwargs, einfo):
        data = {
            "exc": str(exc),
            "task_id": task_id,
            "args": args,
            "kwargs": kwargs,
            **get_base_log_data(),
        }
        project_id = kwargs.get("project_id")
        if project_id:
            project = SampleAnnotationProject.objects.get(id=project_id)
            project.close_process()
            log_to_file(project_id, data)
        log_error.s(data).apply_async(queue=log_queue)


metwork_task = shared_task(base=MetWorkTask)


def log_to_file(project_id, data):
    data = {**data, **get_base_log_data()}
    log_project.s(project_id, data).apply_async(queue=log_queue)


def get_base_log_data():
    return {"time": datetime.now(), "hostname": socket.gethostname()}


@metwork_task
def start_run(project_id):

    project = SampleAnnotationProject.objects.get(id=project_id)
    cache.set(project.process_count_key(), 0, None)
    if project.depth_total > 0:
        project.status_code = SampleAnnotationProject.status.RUNNING
        project.save()
        cache.set(
            "project_ms1_not_init_" + str(project_id), project.ms1_not_init(), None
        )
        molecules_all = [m.id for m in project.molecules.all()]
        cache.set("project_molecules_all_" + str(project_id), molecules_all, None)
        project.add_process(len(molecules_all))
        tsk = group(
            run_reactions_molecule.s(
                molecule_id=m_id,
                project_id=project_id,
                depth_total=0,
                depth_last_match=0,
            )
            for m_id in molecules_all
        )
        tsk()
    else:
        project.status_code = SampleAnnotationProject.status.DONE
        project.save()


@metwork_task
def run_reactions_molecule(molecule_id, project_id, depth_total, depth_last_match):

    project = SampleAnnotationProject.objects.get(id=project_id)
    if project.is_stopped():
        project.close_process()
        return "project stopped"

    grp_tsk = group(
        run_reaction_molecule.s(
            [molecule_id], r.id, project_id, depth_total, depth_last_match
        )
        for r in project.reactions.all()
    )

    project.add_process(project.reactions.count())
    grp_tsk()

    project.close_process()
    return 0


@metwork_task
def run_reaction_molecule(
    reactants_id, reaction_id, project_id, depth_total, depth_last_match
):
    from metabolization.models import Reaction, ReactProcess
    from django.db.models import Count

    project = SampleAnnotationProject.objects.get(id=project_id)
    if project.is_stopped():
        project.close_process()
        return "project stopped"
    reaction = model_from_cache(Reaction, reaction_id)

    log_message = {
        "task": "run_reaction_molecule",
        "opts": {
            "reaction_id": reaction_id,
            "smarts": reaction.smarts,
            "reaction_name": reaction.name,
        },
    }

    reactants_count = len(reactants_id)
    reactants_uniq_count = len(set(reactants_id))

    if reactants_count == reaction.reactants_number:
        reactants = Molecule.objects.filter(id__in=reactants_id)
        log_message["opts"].update(
            {"reactants": {r.id: r.smiles() for r in reactants.all()}}
        )

        rp_search = ReactProcess.objects.annotate(
            reactants_count=Count("reactants")
        ).filter(reaction=reaction, reactants_count=reactants_uniq_count)
        for m in reactants:
            rp_search = rp_search.filter(reactants__id=m.id)

        if rp_search.count() > 0:
            log_message.update({"message": "reaction already processed"})
            rp = rp_search.first()
            rp.wait_run_end()
        else:
            rp = ReactProcess.objects.create(reaction=reaction)
            for m_id in reactants_id:
                rp.reactants.add(Molecule.objects.get(id=m_id))
                # reactants_smiles =
            rp.save()
            rp.run_reaction().refresh_from_db()
        try:
            project.react_processes.add(rp)
            project.save()
        except:
            pass

        log_message["opts"].update(
            {"products": {r.id: r.smiles() for r in rp.products.all()}}
        )
        log_to_file(project_id, log_message)

        mols_ids = [m.id for m in rp.products.all()]

        grp_tsk = group(
            load_molecule.s(m_id, project_id, depth_total, depth_last_match)
            for m_id in mols_ids
        )

        project.add_process(len(mols_ids))
        grp_tsk()

    # if missing 2nd reactant, react with each molecule that has already match
    elif reactants_count == 1 and reaction.reactants_number == 2:
        mols_ids = [m.id for m in project.molecules_init_and_matching()]
        grp_tsk = group(
            run_reaction_molecule.s(
                [reactants_id[0], m_id],
                reaction.id,
                project_id,
                depth_total,
                depth_last_match,
            )
            for m_id in mols_ids
        )

        project.add_process(len(mols_ids))
        grp_tsk()

    project.close_process()
    return 0


@metwork_task
def load_molecule(molecule_id, project_id, depth_total, depth_last_match):

    project = SampleAnnotationProject.objects.get(id=project_id)
    if project.is_stopped():
        project.close_process()
        return "project stopped"

    molecules_ids = cache.get("project_molecules_all_" + str(project_id))

    if not molecule_id in molecules_ids:
        molecules_ids.append(molecule_id)
        molecule = model_from_cache(Molecule, molecule_id)

        try:
            project.molecules.add(molecule)
            project.save()
            added = True
            cache.set("project_molecules_all_" + str(project_id), molecules_ids, None)
        except:
            added = False

        if added:
            project.add_process()
            tsk = evaluate_molecule.apply_async(
                args=[molecule_id, project_id, depth_total, depth_last_match],
                queue=settings.CELERY_RUN_QUEUE,
            )

    project.close_process()
    return 0


@metwork_task
def evaluate_molecule(molecule_id, project_id, depth_total, depth_last_match):
    from fragmentation.models import FragMolSample, FragAnnotationDB

    molecule = model_from_cache(Molecule, molecule_id)
    project = SampleAnnotationProject.objects.get(id=project_id)
    if project.is_stopped():
        project.close_process()
        return "project stopped"

    log_message = {
        "task": "evaluate_molecule",
        "opts": {"molecule_id": molecule_id, "smiles": molecule.smiles(),},
    }

    # check if sample exist with the same mass

    fm_search_ids = []
    ms1 = cache.get("project_ms1_not_init_" + str(project_id))
    am = AdductManager(ion_charge=project.frag_sample.ion_charge)
    adducts_mass = np.array(am.adducts.mass).reshape((1, -1))

    mass_exact = molecule.mass_exact()  # + settings.PROTON_MASS

    diff_mass = abs(ms1[1].reshape((-1, 1)) - mass_exact - adducts_mass)
    mass_tol = mass_exact * float(project.frag_compare_conf.ppm_tolerance) * 10 ** -6
    test = np.where(diff_mass <= mass_tol)
    # fm_search_ids = [(mol_id, adduct), ...]
    fm_search_ids = [
        (ms1[0][pos_idx], am.adducts.index[adduct_idx])
        for pos_idx, adduct_idx in zip(test[0], test[1])
    ]

    adducts = [am.adducts.index[adduct_idx] for adduct_idx in set(test[1])]

    log_message["opts"].update({"fm_search_ids": fm_search_ids, "adducts": adducts})
    log_to_file(project_id, log_message)

    project.add_process()
    tsk = evaluate_molecule_2.apply_async(
        args=[
            molecule.id,
            fm_search_ids,
            adducts,
            project_id,
            depth_total,
            depth_last_match,
        ],
        queue=settings.CELERY_RUN_QUEUE,
    )

    project.close_process()
    return 0


@metwork_task
def evaluate_molecule_2(
    molecule_id, fm_search_ids, adducts, project_id, depth_total, depth_last_match
):
    from fragmentation.models import FragMolSample, FragAnnotationCompare
    from fragmentation.modules import FragSim, FragCompare

    molecule = model_from_cache(Molecule, molecule_id)
    project = SampleAnnotationProject.objects.get(id=project_id)
    if project.is_stopped():
        project.close_process()
        return "project stopped"

    fmsamples = [
        (FragMolSample.objects.get(id=fmsample_id), adduct)
        for fmsample_id, adduct in fm_search_ids
    ]
    ion_ids = [(str(fms.ion_id), adduct) for fms, adduct in fmsamples]

    log_message = {
        "task": "evaluate_molecule_2",
        "opts": {
            "molecule_id": molecule_id,
            "smiles": molecule.smiles(),
            "fm_search_ids": fm_search_ids,
            "ion_ids": ion_ids,
            "adducts": adducts,
        },
    }

    # If mass match, frag molecule and check if frag match

    match = False
    if len(fm_search_ids) > 0:
        # fm_search_ids can be = 0 for ghost metabolites

        ### Fragmentation
        fsim = FragSim(project.frag_sim_conf)
        fmsim = {}
        for adduct in adducts:
            fmsim[adduct] = fsim.frag_molecule(
                molecule, adduct, ion_charge=project.frag_sample.ion_charge
            )

        ### Compare with each sample of same mass
        fcomp = FragCompare(project.frag_compare_conf)
        for fmsample_id, adduct in fm_search_ids:
            fmsample = FragMolSample.objects.get(id=fmsample_id)
            fac_search = FragAnnotationCompare.objects.filter(
                project=project, frag_mol_sample=fmsample, molecule=molecule
            )
            if fac_search.count() == 0:
                fac = FragAnnotationCompare.objects.create(
                    project=project, frag_mol_sample=fmsample, molecule=molecule
                )
                fmcomp = fcomp.compare_frag_mols([fmsim[adduct], fmsample])
                fac.frag_mol_compare = fmcomp
                fac.save()
                match = fmcomp.match

                log_message["opts"].update(
                    {"match": str(match), "cosine": fmcomp.cosine}
                )

    log_to_file(project_id, log_message)
    # Metabolize molecule under conditions

    if ((depth_total + 1) < project.depth_total) and (
        match or (depth_last_match < project.depth_last_match)
    ):

        depth_total += 1
        if match:
            depth_last_match = 0
        else:
            depth_last_match += 1

        project.add_process()
        tsk = run_reactions_molecule.apply_async(
            args=[molecule_id, project_id, depth_total, depth_last_match],
            queue=settings.CELERY_RUN_QUEUE,
        )

    project.close_process()
    return 0


@shared_task
def get_metabolization_network(project_id):
    project = SampleAnnotationProject.objects.get(id=project_id)
    project.get_metabolization_network(task=False, force=True)
    return project_id
