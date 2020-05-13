from .utils import get_env

# Send mail
DEFAULT_FROM_EMAIL = "metwork@pharmacie.parisdescartes.fr"
SERVER_EMAIL = "metwork.dev@gmail.com"
EMAIL_USE_TLS = True
EMAIL_HOST = "smtp.gmail.com"
EMAIL_PORT = 587
EMAIL_HOST_USER = "metwork.dev@gmail.com"
EMAIL_HOST_PASSWORD = get_env("METWORK_EMAIL_HOST_PASSWORD")

GUEST_USER_EMAIL = "metwork.dev@gmail.com"
