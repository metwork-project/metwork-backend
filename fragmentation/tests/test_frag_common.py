# -*- coding: utf-8 -*-
from __future__ import unicode_literals


class FragCommonTests(object):
    @classmethod
    def new_frag_sim(self, name="test"):
        from fragmentation.models import FragSimConf
        from fragmentation.modules import FragSim

        fsc = FragSimConf.objects.create(
            param_path="param/param_output0.log", conf_path="conf/param_config.txt",
        )

        return FragSim(fsc)

    test_data_6 = {
        "smiles": "N=C([NH3+])",
        "result": {
            0: [(18.033826, 0.397148), (28.018175, 0.993644), (45.044725, 98.609208)],
            1: [(18.033826, 4.300341), (28.018175, 1.870827), (45.044725, 93.828831)],
            2: [(18.033826, 4.580884), (28.018175, 6.700656), (45.044725, 88.718460)],
        },
    }

    test_data = {
        "smiles": "N=C([NH3+])",
        "result": {
            0: [
                (18.03382555, 0.3971478062),
                (28.01817548, 0.9936442241),
                (45.04472458, 98.6092079700),
            ],
            1: [
                (18.03382555, 4.30034122),
                (28.01817548, 1.870827298),
                (45.04472458, 93.828831480),
            ],
            2: [
                (18.03382555, 4.580883545),
                (28.01817548, 6.700656456),
                (45.04472458, 88.718460000),
            ],
        },
    }
