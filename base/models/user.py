#https://simpleisbetterthancomplex.com/tutorial/2016/07/22/how-to-extend-django-user-model.html#abstractbaseuser

from __future__ import unicode_literals

from django.db import models
from django.core.mail import send_mail
from django.contrib.auth.models import PermissionsMixin
from django.contrib.auth.base_user import AbstractBaseUser

from .managers import UserManager

class User(AbstractBaseUser, PermissionsMixin):
    email = models.EmailField(
                unique=True)
    username = models.CharField(
                max_length=30,
                blank=True)
    first_name = models.CharField(
                max_length=30,
                blank=True)
    last_name = models.CharField(
                max_length=30,
                blank=True)
    organization = models.CharField(
                max_length=30)
    date_joined = models.DateTimeField(
                auto_now_add=True)
    is_active = models.BooleanField(
                default=True)
    is_staff = models.BooleanField(
                default=False)

    objects = UserManager()

    USERNAME_FIELD = 'email'
    REQUIRED_FIELDS = []

    class Meta:
        verbose_name = 'user'
        verbose_name_plural = 'users'

    class JSONAPIMeta:
        resource_name = "users"

    def __str__(self):
        return 'user' #self.email

    def get_short_name(self):
        return self.email

    def email_user(self, subject, message, from_email="metwork@pharmacie.parisdescartes.fr", **kwargs):
        '''
        Sends an email to this User.
        '''
        send_mail(subject, message, from_email, [self.email], **kwargs)
