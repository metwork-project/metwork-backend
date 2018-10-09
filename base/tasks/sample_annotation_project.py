from __future__ import absolute_import, unicode_literals
from celery import shared_task, group, chain
from base.models import SampleAnnotationProject, Molecule
from base.modules.cache_management import model_from_cache
from django.core.cache import cache
from django.conf import settings
import numpy as np

@shared_task
def start_run(project_id):
    p = SampleAnnotationProject.objects.get(id = project_id)
    cache.set( p.process_count_key(), 0, None )
    if p.depth_total > 0:
        p.status_code = SampleAnnotationProject.status.RUNNING
        p.save()
        cache.set( 'project_ms1_not_init_' + str(project_id), p.ms1_not_init(), None )
        molecules_all = [ m.id for m in p.molecules.all() ]
        cache.set('project_molecules_all_' + str(project_id), molecules_all, None)
        p.add_process(len(molecules_all))
        tsk = group(\
                run_reactions_molecule.s( m_id, project_id, 0, 0 )\
                for m_id in molecules_all )
        tsk()
    else:
        p.status_code = SampleAnnotationProject.status.DONE
        p.save()

@shared_task
def run_reactions_molecule(molecule_id, project_id, depth_total, depth_last_match ):

    m = model_from_cache( Molecule, molecule_id )
    p = model_from_cache( SampleAnnotationProject, project_id )

    grp_tsk = group(\
                run_reaction_molecule.s([molecule_id], r.id, project_id, depth_total, depth_last_match) \
                for r in p.reactions_conf.reactions.all() )

    p.add_process(p.reactions_conf.reactions.count())
    grp_tsk()

    p.close_process()
    return 0

@shared_task
def run_reaction_molecule(reactants_id, reaction_id, project_id, depth_total, depth_last_match ):
    from metabolization.models import Reaction, ReactProcess
    from django.db.models import Count

    p = model_from_cache( SampleAnnotationProject, project_id )
    r = model_from_cache( Reaction,reaction_id )

    reactants_count = len(reactants_id)
    reactants_uniq_count = len(set(reactants_id))

    if reactants_count == r.reactants_number:
        reactants = Molecule.objects\
                        .filter(id__in = reactants_id)

        rp_search = ReactProcess.objects\
                        .annotate(reactants_count=Count('reactants'))\
                        .filter(
                            reaction = r,
                            reactants_count = reactants_uniq_count,
                            method = r.method_to_apply())
        for m in reactants:
            rp_search = rp_search.filter(reactants__id = m.id)

        if rp_search.count() > 0:
            rp = rp_search.first()
            rp.wait_run_end()
        else:
            rp = ReactProcess.objects.create(reaction = r)
            for m_id in reactants_id:
                rp.reactants.add(Molecule.objects.get( id = m_id ))
            rp.save()
            rp.run_reaction().refresh_from_db()
        try:
            p.react_processes.add(rp)
            p.save()
        except:
            pass

        mols_ids = [m.id for m in rp.products.all() ]

        grp_tsk = group( load_molecule.s(m_id, project_id, depth_total, depth_last_match) for m_id in mols_ids )

        p.add_process(len(mols_ids))
        grp_tsk()

    # if missing 2nd reactant, react with each molecule that has already match
    elif reactants_count == 1 and r.reactants_number == 2:
        mols_ids = [m.id for m in p.molecules_init_and_matching() ]
        grp_tsk = group(\
            run_reaction_molecule.s(\
                [reactants_id [0], m_id], r.id, project_id, depth_total, depth_last_match ) \
                for m_id in mols_ids )

        p.add_process( len(mols_ids) )
        grp_tsk()

    p.close_process()
    return 0

@shared_task
def load_molecule(molecule_id, project_id, depth_total, depth_last_match):

    p = model_from_cache( SampleAnnotationProject, project_id )
    molecules_ids = cache.get('project_molecules_all_' + str(project_id))

    if not molecule_id in molecules_ids :
        molecules_ids.append(molecule_id)
        m = model_from_cache( Molecule, molecule_id )

        try:
            p.molecules.add(m)
            p.save()
            added = True
            cache.set('project_molecules_all_' + str(p.id), molecules_ids, None)
        except:
            added = False

        if added:
            p.add_process()
            tsk = evaluate_molecule.apply_async( args= [molecule_id, project_id, depth_total, depth_last_match], queue = settings.CELERY_RUN_QUEUE)

    p.close_process()
    return 0

@shared_task
def evaluate_molecule(molecule_id, project_id, depth_total, depth_last_match):
    from fragmentation.models import FragMolSample, FragAnnotationDB

    m = model_from_cache( Molecule, molecule_id )
    p = model_from_cache( SampleAnnotationProject, project_id )

    # check if sample exist with the same mass

    project_frag_sample = cache.get_or_set(
        'project_frag_sample_' + str(p.id),
        p.frag_sample )

    fm_search_ids = []
    #for mass_exact in m.mass_exact_isotopes():
        #mass_exact += 1.007
    #for mass_exact in [m.mass_exact()]:
    mass_exact = m.mass_exact() + settings.PROTON_MASS
    mass_var = float(p.frag_compare_conf.ppm_tolerance) * 10**-6
    mass_window = (
        mass_exact * ( 1 - mass_var) ,
        mass_exact * (1 + mass_var) )

    ms1 = cache.get( 'project_ms1_not_init_' + str(project_id) )
    pos_id_min, pos_id_max = ( np.searchsorted(ms1[1], mw) for mw in mass_window )
    fm_search_ids = ms1[0][pos_id_min:pos_id_max]

    # ===> Tautomerization not used
    if False: #fm_search.count() > 0:
        maj_taut = m.major_tautomers()
        maj_taut_to_process = maj_taut
        for mt in maj_taut:
            if not mt in p.molecules.all():
                try:
                    p.molecules.add(mt)
                    p.save()
                    maj_taut_to_process.append(mt)
                except:
                    pass
        mol_ids = list( {m.id} | {mol.id for mol in maj_taut_to_process} )
    #else:
        mol_ids = [m.id]
        fm_search_ids += [fms.id for fms in fm_search]
    # =====

    p.add_process()
    tsk = evaluate_molecule_2.apply_async( args= [m.id, fm_search_ids, project_id, depth_total, depth_last_match], queue = settings.CELERY_RUN_QUEUE)

    p.close_process()
    return 0

@shared_task
def evaluate_molecule_2(molecule_id, fm_search_ids, project_id, depth_total, depth_last_match):
    from fragmentation.models import FragMolSample, FragAnnotationCompare
    from fragmentation.modules import FragSim, FragCompare

    m = model_from_cache( Molecule, molecule_id )
    p = model_from_cache( SampleAnnotationProject, project_id )

    # If mass match, frag molecule and check if frag match

    match = False
    if len(fm_search_ids) > 0:
    # fm_search_ids can be = 0 for ghost metabolites

        ### Fragmentation
        fsim = FragSim(p.frag_sim_conf)
        fmsim = fsim.frag_molecule(m)

        ### Compare with each sample of same mass
        fcomp = FragCompare(p.frag_compare_conf)
        for fmsample_id in fm_search_ids:
            fmsample = FragMolSample.objects.get(id = fmsample_id)
            fac_search = FragAnnotationCompare.objects.filter(
                project = p,
                frag_mol_sample = fmsample,
                molecule = m)
            if fac_search.count() == 0:
                fac = FragAnnotationCompare.objects.create(
                    project = p,
                    frag_mol_sample = fmsample,
                    molecule = m)
                fmcomp = fcomp.compare_frag_mols([fmsim, fmsample])
                fac.frag_mol_compare = fmcomp
                fac.save()
                match = fmcomp.match

    # Metabolize molecule under conditions

    if ( (depth_total + 1) < p.depth_total ) \
        and ( match or (depth_last_match < p.depth_last_match) ):

        depth_total += 1
        if match:
            depth_last_match = 0
        else:
            depth_last_match += 1

        p.add_process()
        tsk = run_reactions_molecule.apply_async( args= [molecule_id, project_id, depth_total, depth_last_match], queue = settings.CELERY_RUN_QUEUE)

    p.close_process()
    return 0
