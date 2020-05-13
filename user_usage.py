from base.models import User

print "name, reactions, frag_samples, projects"

print [
    (
        u.username,
        u.reaction_set.count(),
        u.fragsample_set.count(),
        u.project_set.count(),
    )
    for u in User.objects.all()
]
