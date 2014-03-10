from test_main import compare_dicts
from interfits.interfits import InterFits, PrintLog
from interfits.ledafits import LedaFits
import os
import numpy as np
import pylab as plt
import time

pp = PrintLog()
h1 = pp.h1
h2 = pp.h2

import pprint
ppp = pprint.PrettyPrinter(indent=4)
ppp = ppp.pprint

def test_dada():

    ddd = LedaFits('data/test.dada')

    h2("HEADER: COMMON")
    ppp(ddd.h_common)

    h2("HEADER: PARAMS")
    ppp(ddd.h_params)

    h2("HEADER: UV DATA")
    ppp(ddd.h_uv_data)

    h2("HEADER: ARRAY GEOMETRY")
    ppp(ddd.h_array_geometry)

    h2("HEADER: ANTENNA")
    ppp(ddd.h_antenna)

    h2("HEADER: FREQUENCY")
    ppp(ddd.h_frequency)

    h2("HEADER: SOURCE")
    ppp(ddd.h_frequency)

    h2("DATA: FREQUENCY")
    ppp(ddd.d_frequency)


if __name__ == '__main__':

    test_dada()