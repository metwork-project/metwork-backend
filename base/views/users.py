from __future__ import unicode_literals

from base.views.model_auth import ModelAuthViewSet
#from django.contrib.auth.models import User
#from base.models import User
from rest_framework import serializers
from django.contrib.auth import get_user_model
from django.http import HttpResponse, JsonResponse
from django.views.decorators.http import require_http_methods
from django.views.decorators.csrf import csrf_exempt
from base.forms import RegistrationForm, ResetPasswordForm, BetaTestRegistrationForm
import json

class UserSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = get_user_model()
        fields = ('email', 'username', 'organization', )

class UserViewSet(ModelAuthViewSet):
    queryset = get_user_model().objects.all()
    serializer_class = UserSerializer

@require_http_methods(["POST"])
@csrf_exempt
def betatest_register(request):
    """
    API endpoint to register a new user.
    """
    try:
        payload = json.loads(request.body)
    except ValueError:
        return JsonResponse({"error": "Unable to parse request body"}, status=400)

    form = BetaTestRegistrationForm(payload)

    if form.is_valid():
        get_user_model().betatest_register(
            {'email': form.cleaned_data["email"],
            'name': form.cleaned_data["name"],
            'organization': form.cleaned_data["organization"]})
        return JsonResponse({"success": "User registered for betatest."}, status=201)

    return HttpResponse(form.errors.as_json(), status=400, content_type="application/json")

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

        #try:
        user.email_user(
            subject = 'Metwork account created',
            message = 'Welcome to Metwork {0}!\n\nYour account has been created.'.format(username))
        #except:
        #    pass

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

    form = ResetPasswordForm(payload)

    if form.is_valid():
        user = get_user_model()\
                .objects.get(email = form.cleaned_data["email"])

        temp_pwd = get_user_model().objects.make_random_password(length=10)

        user.set_password(temp_pwd)
        user.save()

        user.email_user(
            subject = 'MetWork beta test password',
            message = 'Welcome to MetWork {0} !\n\nYou can now begin beta test on https://metwork.pharmacie.parisdescartes.fr\n\nYour email for login is : {1}\nYour password : {2}\n\nHave fun !\nThe MetWork Team'.format(user.first_name, user.email, temp_pwd))

        return JsonResponse({"success": "Reset paswword email sent"}, status=201)

    return HttpResponse(form.errors.as_json(), status=400, content_type="application/json")
