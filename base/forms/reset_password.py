# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django import forms
#from base.models import User
from django.contrib.auth import get_user_model

class ResetPasswordForm(forms.Form):

    email = forms.EmailField()
 
    def clean_email(self):
        email = self.cleaned_data["email"]
 
        try:
            user = get_user_model().objects.get(email=email)
            return user
        except get_user_model().DoesNotExist:
            raise forms.ValidationError(
                "The email %s does not exist." % email)

