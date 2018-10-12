"""metwork URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.11/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.conf.urls import url, include
from django.contrib import admin
from base.views import APIStatusViewSet, UserViewSet, \
    MoleculeViewSet, ProjectViewSet, \
    register, password_reset #, get_status
from metabolization.views import ReactionViewSet
from fragmentation.views import FragSampleViewSet, FragAnnotationViewSet, FragCompareConfViewSet
from rest_framework import routers, serializers, viewsets
from rest_framework.authtoken.views import obtain_auth_token
from django.conf import settings
from django.conf.urls.static import static

# Routers provide an easy way of automatically determining the URL conf.
router = routers.DefaultRouter(trailing_slash=False)
router.register(r'api-statuses', APIStatusViewSet)
router.register(r'users', UserViewSet)
router.register(r'molecules', MoleculeViewSet)
router.register(r'projects', ProjectViewSet)
router.register(r'fragsamples', FragSampleViewSet)
router.register(r'frag-annotations', FragAnnotationViewSet)
router.register(r'reactions', ReactionViewSet)
router.register(r'frag-compare-confs', FragCompareConfViewSet)

# Wire up our API using automatic URL routing.
# Additionally, we include login URLs for the browsable API.

urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'^', include(router.urls)),
    url(r'^api-auth-token/', obtain_auth_token),
    url(r'^api-auth/', include('rest_framework.urls', namespace='')),
    url(r'^api-register/', register),
    url(r'^api-password-reset/', password_reset),
    # url(r'^api-status/', get_status),
]
