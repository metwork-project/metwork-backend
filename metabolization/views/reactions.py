# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.db.models import Q, Prefetch, Count, Min
from base.views.model_auth import ModelAuthViewSet, IsOwnerOrPublic
from collections import defaultdict
from metabolization.models import Reaction
from rest_framework import serializers
from rest_framework.decorators import action
from rest_framework.response import Response
from django.contrib.auth import get_user_model
from rest_framework.parsers import JSONParser
from django.http import JsonResponse
from base.models import Molecule
from base.modules import JSONSerializerField, ChemDoodle, TagViewMethods
from base.modules.queryset import FilteredQueryset


class ReactionSerializer(serializers.ModelSerializer):
    user = serializers.PrimaryKeyRelatedField(
        queryset=get_user_model().objects.all()
        # ,allow_null = True
    )

    class Meta:
        model = Reaction
        fields = (
            "name",
            "description",
            "tags_list",
            "user",
            "user_id",
            "user_name",
            "reactants_number",
            "has_no_project",
            "status_code",
            "smarts",
            "chemdoodle_json",
            "chemdoodle_json_error",
        )

    chemdoodle_json = JSONSerializerField()


class ReactionQueryset(FilteredQueryset):

    def get_all(self, queryset_):
        return queryset_.all_reactions()

    def get_notselected(self, queryset_):
        return queryset_.reactions_not_selected()

    def filter_other(self):
        queryset = self.queryset

        for key, value in self.request.query_params.lists():

            if key == "filter[text]":
                text = value[0]
                if text:
                    queryset = queryset.filter(
                        Q(name__icontains=text) | Q(tags__name__icontains=text)
                    )
            elif key == "filter[my]":
                my = value[0].lower() == "true"
                if my:
                    queryset = queryset.filter(user=self.request.user)
            elif key == "filter[user]" and not my:
                user_text = value[0]
                if user_text:
                    queryset = queryset.filter(user__username__icontains=user_text)

        self.queryset = queryset


class ReactionViewSet(ModelAuthViewSet, TagViewMethods):
    queryset = Reaction.objects.all().order_by("name")
    serializer_class = ReactionSerializer
    permission_classes = (IsOwnerOrPublic,)

    def get_queryset(self):
        queryset = ReactionQueryset(self.request, Reaction).queryset
        return queryset.order_by("name")

    def list(self, request, *args, **kwargs):
        response = super().list(request, *args, **kwargs)
        queryset = self.filter_queryset(self.get_queryset())
        ids = [reaction.id for reaction in queryset.all()]
        response.data["meta"]["ids"] = ids
        return response

    def create(self, request, *args, **kwargs):
        request.data["user"] = request.user.id
        return super().create(request, *args, **kwargs)

    def update(self, request, *args, **kwargs):
        if "user" in request.data:
            del request.data["user"]
        return super().update(request, *args, **kwargs)

    @action(detail=False, methods=["post"])
    def uploadfile(self, request):
        req_data = request.data
        try:
            fs = Reaction.import_file(
                file_object=req_data["file_data"],
                user=request.user,
                name=req_data["name"],
                description=req_data["description"],
            )
            return Response({"status": "ok"})
        except:
            return Response({"error": "unkown error while upploading file"})

    @action(detail=True, methods=["get"])
    def get_image(self, request, pk=None):
        reaction = self.get_object()
        return Response({"image": reaction.get_image()})

    @action(detail=True, methods=["patch"])
    def load_smarts(self, request, pk=None):
        try:
            reaction = self.get_object()
            serializer = self.serializer_class(reaction)
            data = JSONParser().parse(request)
            reaction.load_smarts(data["smarts"])
            return JsonResponse({"success": reaction.chemdoodle_json})
        except:
            return JsonResponse({"error": "error import smarts"})

    @action(detail=True, methods=["post"])
    def run_reaction(self, request, pk=None):
        data = JSONParser().parse(request)
        chemdoodle_json = data["reactants"]["chemdoodle_json"]
        if "m" in chemdoodle_json:
            cd = ChemDoodle()
            reactants = [cd.json_to_mol(mol_json) for mol_json in chemdoodle_json["m"]]
            reactants_smiles = reactants[0].smiles()
            if "." in reactants_smiles:
                reactants = [
                    Molecule.load_from_smiles(sm) for sm in reactants_smiles.split(".")
                ]
            r = self.get_object()
            rp = r.run_reaction(reactants)
            response = {
                "reactants": [r.chemdoodle_json for r in rp.reactants.all()],
                "products": [p.chemdoodle_json for p in rp.products.all()],
            }
            return JsonResponse(response)
        else:
            return JsonResponse({"error": "test"})
