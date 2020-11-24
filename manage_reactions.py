from metabolization.models import Reaction
from base.models import User

FILE_PATH = "reactions_export.tsv"
SEPARATOR = "\t"


def export_reactions(file_path=FILE_PATH):
    with open(file_path, "w") as fw:
        for r in Reaction.objects.filter(status_code=Reaction.status.ACTIVE):
            line = SEPARATOR.join([r.name, r.smarts]) + "\n"
            fw.write(line)


def import_reactions(file_path=FILE_PATH, email="metwork.dev@gmail.com"):
    user = User.objects.get(email=email)
    with open(file_path, "r") as fr:
        lines = r = fr.readlines()
        for line in lines:
            try:
                r_raw = line.replace("\n", "").split(SEPARATOR)
                r = Reaction(name=r_raw[0], user=user)
                r.save()
                r.load_smarts(r_raw[1])
                r.status_code = Reaction.status.ACTIVE
                r.save()
            except:
                print("error importing {}".format(r_raw[0]))
