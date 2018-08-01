# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.test import TransactionTestCase
from fragmentation.modules import FragCompare
from fragmentation.models import *
from django.core.exceptions import FieldError

class FragMolMolCompare(TransactionTestCase):

    def test_2_fragmols(self):
        fconf = FragCompareConf.objects.create()
        fcomp = FragCompare(fconf)
        fm1 = FragMol.objects.create()
        fm2 = FragMol.objects.create()
        fm3 = FragMol.objects.create()
        fmc = FragMolCompare.objects.create()
        with self.assertRaises(FieldError): 
            fmc.frag_mols.add(fm1)
            fmc.save()   

        try:
            fmc.frag_mols.add(fm2)
            fmc.save()  
        except:
            self.fail("fmc.save() fail")

        try:
            fcomp.compare_frag_mols([fm1, fm2])
        except:
            self.fail("fcomp.compare_frag_mols([fm, fm]) fail")
       
        with self.assertRaises(FieldError): 
            fcomp.compare_frag_mols([fm1])
        with self.assertRaises(FieldError): 
            fcomp.compare_frag_mols([fm1, fm2, fm3])

