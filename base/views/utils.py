from base.views.model_auth import ModelAuthViewSet

class MetaIdsViewSet(ModelAuthViewSet):

    def list(self, request, *args, **kwargs):
        response = super().list(request, *args, **kwargs)
        queryset = self.filter_queryset(self.get_queryset())
        ids = [item.id for item in queryset.all()]
        response.data["meta"]["ids"] = ids
        return response