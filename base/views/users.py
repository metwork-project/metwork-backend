from __future__ import unicode_literals

from base.views.model_auth import ModelAuthViewSet
#from django.contrib.auth.models import User
#from base.models import User
from rest_framework import serializers
from rest_framework.decorators import list_route, detail_route
from django.contrib.auth import get_user_model
from django.http import HttpResponse, JsonResponse
from rest_framework.response import Response
from django.views.decorators.http import require_http_methods
from django.views.decorators.csrf import csrf_exempt
from base.forms import RegistrationForm, ResetPasswordForm, BetaTestRegistrationForm
import json
from rest_framework.parsers import JSONParser
from django.conf import settings

class UserSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = get_user_model()
        fields = ('email', 'username', 'organization', )

class UserViewSet(ModelAuthViewSet):
    queryset = get_user_model().objects.all()
    serializer_class = UserSerializer

    def get_queryset(self):
        return get_user_model().objects.filter(id = self.request.user.id)

    def update(self, *args, **kwargs):
        data = args[0].data
        if  data['email'] != settings.GUEST_USER_EMAIL:
            super(UserViewSet, self).update( *args, **kwargs )
        return Response( UserSerializer( get_user_model().objects.get(id=data['id']) ).data )

    @detail_route(methods=['patch'])
    def change_password(self, request, pk=None):
        if  self.request.user.email != settings.GUEST_USER_EMAIL:
            user = self.get_object()
            password = JSONParser().parse(request)['password']
            user.set_password(password)
            user.save()
            return Response({'success': 'Password changed'})

@require_http_methods(["POST"])
@csrf_exempt
def register(request):
    """
    API endpoint to register a new user.
    """
    try:
        payload = json.loads(request.body)
    except ValueError:
        return JsonResponse({"error": "Unable to parse request body"}, status=400)

    form = RegistrationForm(payload)

    if form.is_valid():
        username = form.cleaned_data["username"]
        user = get_user_model().objects.create_user(
                email = form.cleaned_data["email"],
                username = username,
                organization = form.cleaned_data["organization"],
                password = form.cleaned_data["password"])
        user.save()

        try:
            user.email_user(
                subject = 'Metwork account created',
                message = 'Welcome to Metwork {0}!\n\nYour account has been created.'.format(username))
        except:
            pass

        return JsonResponse({"success": "User registered."}, status=201)

    return HttpResponse(form.errors.as_json(), status=400, content_type="application/json")

@require_http_methods(["POST"])
@csrf_exempt
def password_reset(request):

    """
    In case of forgot password
    """
    try:
        payload = json.loads(request.body)
    except ValueError:
        return JsonResponse({"error": "Unable to parse request body"}, status=400)

    print('payload',payload)

    form = ResetPasswordForm(payload)

    if form.is_valid():
        print('form',form.cleaned_data)
        print("email", form.cleaned_data["email"])
        user = get_user_model()\
                .objects.get(email = form.cleaned_data["email"])

        temp_pwd = get_user_model().objects.make_random_password(length=10)

        user.set_password(temp_pwd)
        user.save()

        user.email_user(
            subject = '[MetWork] Password reset',
            message = 'Your new password for MetWork is : {0}\n\nThe MetWork Team'.format(temp_pwd) )

        return JsonResponse({"success": "Reset password email sent"}, status=201)

    return HttpResponse(form.errors.as_json(), status=400, content_type="application/json")
