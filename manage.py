#!/usr/bin/env python
import os
import sys

if __name__ == "__main__":
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "metwork_backend.settings")
    try:
        from django.core.management import execute_from_command_line
    except ImportError:
        # The above import may fail for some other reason. Ensure that the
        # issue is really that Django is missing to avoid masking other
        # exceptions on Python 2.
        try:
            import django
        except ImportError:
            raise ImportError(
                "Couldn't import Django. Are you sure it's installed and "
                "available on your PYTHONPATH environment variable? Did you "
                "forget to activate a virtual environment?"
            )
        raise

    first = not "METWORK_DEV_RUNNING" in os.environ
    if sys.argv[1] == 'runserver' and first:
        os.environ["METWORK_DEV_RUNNING"] = "True"
        print ("""
    __ _       __           __  
   /  |/  /__  / /| |     / /___  _____/ /__
  / /|_/ / _ \/ __/ | /| / / __ \/ ___/ //_/
 / /  / /  __/ /_ | |/ |/ / /_/ / /  / ,<   
/_/  /_/\___/\__/ |__/|__/\____/_/  /_/|_|  

**********      Dev Server    ************

        """)
    execute_from_command_line(sys.argv)
    if "METWORK_DEV_RUNNING" in os.environ:
        os.environ.pop("METWORK_DEV_RUNNING")
