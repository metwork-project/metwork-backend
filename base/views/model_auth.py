from rest_framework import viewsets
from rest_framework.authentication import TokenAuthentication, SessionAuthentication
from rest_framework.permissions import IsAuthenticated, BasePermission
from rest_framework import permissions

class ModelAuthViewSet(viewsets.ModelViewSet):
    authentication_classes = (TokenAuthentication, SessionAuthentication)
    permission_classes = (IsAuthenticated,)

class IsOwnerOrPublic(BasePermission):

    def has_object_permission(self, request, view, obj):

        if view.action in ['run_reaction']:
            return True

        if type(obj.user).__name__ == 'method':
            user = obj.user()
        else:
            user = obj.user  

        if  user == request.user:
            return True

        if request.method in permissions.SAFE_METHODS:
            try:
                return obj.public
            except:
                pass
            try:
                return obj.is_public()
            except:
                pass
            return False
