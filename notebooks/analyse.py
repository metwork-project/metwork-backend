from django.contrib.auth import get_user_model
from base.models import SampleAnnotationProject
from fragmentation.models import FragSample
from metabolization.models import Reaction

def write_json(filename, data):
    import json
    from pathlib import Path
    from django.core.serializers.json import DjangoJSONEncoder

    json_values = json.dumps(data, cls=DjangoJSONEncoder)
    Path("{}.json".format(filename)).write_text(json_values)


write_json("projects", list(SampleAnnotationProject.objects.all().values()))


data = [
    {
        "id": p.id,
        "matching_count": p.molecules_matching_count(),
        "reactions_count": p.all_reactions().count(),
    }
    for p in SampleAnnotationProject.objects.all()
]
write_json("projects_info", data)


write_json("fragsamples", list(FragSample.objects.all().values()))

data = [
    {"id": fs.id, "annotations_count": fs.annotations_count()}
    for fs in FragSample.objects.all()
]
write_json("fragsamples_info", data)

write_json("reactions", list(Reaction.objects.all().values()))

User = get_user_model()
write_json("users", list(User.objects.all().values()))
