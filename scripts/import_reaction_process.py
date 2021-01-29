import json
from pathlib import Path
from metabolization.models import Reaction, ReactProcess


username = "Romain MAGNY"

for reaction in Reaction.objects.filter(status_code=30).filter(user__username=username):
    # print(dir(reaction))

    # query = ReactProcess.objects.filter(reaction__status_code=30).filter(
    #     reaction__user__username=username
    # )
    # query = ReactProcess.objects.filter(reaction__status_code=30).filter(
    #     reaction__user__username=username
    # )
    path = Path("../reaction_analysis/RM/{}.json".format(reaction.id))

    data = {}

    f = path.open("w")

    data["reaction_info"] = {
        "id": reaction.id,
        "name": reaction.name,
        "smarts": reaction.smarts,
        "reactants_number": reaction.reactants_number,
    }
    # f.write(json.dumps(reaction_info) + ',\n')

    # f.write("[\n")

    # first = True

    data["data"] = []

    for rp in reaction.reactprocess_set.all():

        line = {
            "id": rp.id,
            "reactants": [mol.smiles() for mol in rp.reactants.all()],
            "products": [mol.smiles() for mol in rp.products.all()],
        }
        data["data"].append(line)
        # if not first:
        #     f.write(",\n")
        # else:
        #     first = False
        # f.write(json.dumps(line))

    # f.write("\n]")
    path.write_text(json.dumps(data))

    # f.write("{\n")
