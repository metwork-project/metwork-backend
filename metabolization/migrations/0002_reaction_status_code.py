# Generated by Django 2.0.7 on 2018-09-18 08:39

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('metabolization', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='reaction',
            name='status_code',
            field=models.PositiveSmallIntegerField(db_index=True, default=0),
        ),
    ]