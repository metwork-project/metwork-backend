# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django import forms

# from base.models import User
from django.contrib.auth import get_user_model


class RegistrationForm(forms.Form):
    """
    Register a new user.
    """

    email = forms.EmailField()
    username = forms.CharField(max_length=30, min_length=3, label="User name")
    organization = forms.CharField(max_length=30, min_length=3, label="Organization")
    password = forms.CharField(
        widget=forms.PasswordInput(), min_length=5, label="Password"
    )
    confirm_password = forms.CharField(
        widget=forms.PasswordInput(), min_length=5, label="Confirm Password"
    )

    def clean_email(self):
        email = self.cleaned_data["email"]

        try:
            user = get_user_model().objects.get(email=email)
        except get_user_model().DoesNotExist:
            return email

        raise forms.ValidationError("The email %s is already taken." % email)

    def clean(self):
        """
        Make sure that the two passwords match.
        """
        password = self.cleaned_data.get("password", None)
        confirm_password = self.cleaned_data.get("confirm_password", None)

        if password == confirm_password:
            return self.cleaned_data

        raise forms.ValidationError("The passwords do not match.")


class BetaTestRegistrationForm(forms.Form):
    """
    Register a new user.
    """

    email = forms.EmailField()
    name = forms.CharField(max_length=30, min_length=3, label="Name")
    organization = forms.CharField(max_length=30, min_length=3, label="Organization")

    def clean(self):

        return self.cleaned_data
