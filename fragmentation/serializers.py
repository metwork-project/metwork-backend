from rest_framework import serializers
from fragmentation.models import FragSimConf

class FragSimConfSerializer(serializers.ModelSerializer):
    class Meta:
        model = FragSimConf
        fields = ('name','threshold', )

