# -*- coding: utf-8 -*-
# Generated by Django 1.11.10 on 2018-06-25 14:06
from __future__ import unicode_literals

import base.modules.conf_management
from decimal import Decimal
from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        ('contenttypes', '0002_remove_content_type_name'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('base', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='FragAnnotation',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
            ],
            options={
                'manager_inheritance_from_future': True,
            },
        ),
        migrations.CreateModel(
            name='FragCompareConf',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('cosine_threshold', models.DecimalField(decimal_places=3, default=0.18, max_digits=4)),
                ('ppm_tolerance', models.DecimalField(decimal_places=4, default=5, max_digits=6)),
            ],
            bases=(base.modules.conf_management.ConfManagement, models.Model),
        ),
        migrations.CreateModel(
            name='FragMol',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('mass', models.DecimalField(db_index=True, decimal_places=10, default=0, max_digits=16)),
            ],
            options={
                'manager_inheritance_from_future': True,
            },
        ),
        migrations.CreateModel(
            name='FragMolAttribute',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('title', models.CharField(default='', max_length=32)),
                ('value', models.CharField(default='', max_length=32)),
                ('position', models.PositiveSmallIntegerField(default=0)),
            ],
        ),
        migrations.CreateModel(
            name='FragMolCompare',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('match', models.BooleanField(default=False)),
                ('cosine', models.DecimalField(decimal_places=3, default=Decimal('0'), max_digits=4)),
                ('energies', models.CharField(default='', max_length=32)),
                ('num_frag_match', models.IntegerField(default=0)),
                ('frag_compare_conf', models.ForeignKey(blank=True, default=None, null=True, on_delete=django.db.models.deletion.PROTECT, to='fragmentation.FragCompareConf')),
            ],
        ),
        migrations.CreateModel(
            name='FragMolPeak',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('energy', models.SmallIntegerField(default=0)),
                ('mz', models.DecimalField(decimal_places=10, max_digits=26)),
                ('intensity', models.DecimalField(decimal_places=10, max_digits=26)),
            ],
            options={
                'ordering': ('mz',),
            },
        ),
        migrations.CreateModel(
            name='FragSample',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(default='', max_length=128)),
                ('file_name', models.CharField(default='', max_length=255)),
                ('description', models.CharField(blank=True, default='', max_length=255, null=True)),
                ('ions_total', models.PositiveSmallIntegerField(db_index=True, default=0)),
                ('status_code', models.PositiveIntegerField(db_index=True, default=0)),
                ('user', models.ForeignKey(default=None, on_delete=django.db.models.deletion.PROTECT, to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'ordering': ('name',),
            },
        ),
        migrations.CreateModel(
            name='FragSimConf',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('threshold', models.DecimalField(decimal_places=6, default=0.003, max_digits=7)),
                ('param_path', models.CharField(default='param/param_output0.log', max_length=255)),
                ('conf_path', models.CharField(default='conf/param_config.txt', max_length=255)),
            ],
            bases=(base.modules.conf_management.ConfManagement, models.Model),
        ),
        migrations.CreateModel(
            name='FragAnnotationCompare',
            fields=[
                ('fragannotation_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='fragmentation.FragAnnotation')),
            ],
            options={
                'manager_inheritance_from_future': True,
            },
            bases=('fragmentation.fragannotation',),
        ),
        migrations.CreateModel(
            name='FragAnnotationDB',
            fields=[
                ('fragannotation_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='fragmentation.FragAnnotation')),
                ('name', models.CharField(default='', max_length=64)),
                ('db_source', models.CharField(default='unkown', max_length=64)),
                ('db_id', models.CharField(default='', max_length=64)),
            ],
            options={
                'manager_inheritance_from_future': True,
            },
            bases=('fragmentation.fragannotation',),
        ),
        migrations.CreateModel(
            name='FragMolSample',
            fields=[
                ('fragmol_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='fragmentation.FragMol')),
                ('ion_id', models.IntegerField(default=0)),
                ('frag_sample', models.ForeignKey(default=None, null=True, on_delete=django.db.models.deletion.CASCADE, to='fragmentation.FragSample')),
            ],
            options={
                'manager_inheritance_from_future': True,
            },
            bases=('fragmentation.fragmol',),
        ),
        migrations.CreateModel(
            name='FragMolSim',
            fields=[
                ('fragmol_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='fragmentation.FragMol')),
                ('param_hash', models.CharField(default='', max_length=32)),
                ('conf_hash', models.CharField(default='', max_length=32)),
                ('status_code', models.SmallIntegerField(default=0)),
                ('frag_sim_conf', models.ForeignKey(default=None, null=True, on_delete=django.db.models.deletion.CASCADE, to='fragmentation.FragSimConf')),
                ('molecule', models.ForeignKey(default=None, null=True, on_delete=django.db.models.deletion.PROTECT, to='base.Molecule')),
            ],
            options={
                'manager_inheritance_from_future': True,
            },
            bases=('fragmentation.fragmol',),
        ),
        migrations.AddField(
            model_name='fragmolpeak',
            name='frag_mol',
            field=models.ForeignKey(default=None, on_delete=django.db.models.deletion.CASCADE, to='fragmentation.FragMol'),
        ),
        migrations.AddField(
            model_name='fragmolcompare',
            name='frag_mols',
            field=models.ManyToManyField(blank=True, default=None, to='fragmentation.FragMol'),
        ),
        migrations.AddField(
            model_name='fragmolattribute',
            name='frag_mol',
            field=models.ForeignKey(default=None, on_delete=django.db.models.deletion.CASCADE, to='fragmentation.FragMol'),
        ),
        migrations.AddField(
            model_name='fragmol',
            name='polymorphic_ctype',
            field=models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='polymorphic_fragmentation.fragmol_set+', to='contenttypes.ContentType'),
        ),
        migrations.AddField(
            model_name='fragannotation',
            name='molecule',
            field=models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, to='base.Molecule'),
        ),
        migrations.AddField(
            model_name='fragannotation',
            name='polymorphic_ctype',
            field=models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='polymorphic_fragmentation.fragannotation_set+', to='contenttypes.ContentType'),
        ),
        migrations.AddField(
            model_name='fragannotationcompare',
            name='frag_mol_compare',
            field=models.ForeignKey(default=None, null=True, on_delete=django.db.models.deletion.PROTECT, to='fragmentation.FragMolCompare'),
        ),
        migrations.AddField(
            model_name='fragannotationcompare',
            name='project',
            field=models.ForeignKey(default=None, on_delete=django.db.models.deletion.CASCADE, to='base.Project'),
        ),
        migrations.AddField(
            model_name='fragannotation',
            name='frag_mol_sample',
            field=models.ForeignKey(default=None, on_delete=django.db.models.deletion.CASCADE, to='fragmentation.FragMolSample'),
        ),
    ]
