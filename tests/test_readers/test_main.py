import sys
import numpy as np
from interfits.ledafits import *


def compare_dicts(dict_a, dict_b):
    """ Compare two dictionaries to confirm they contain same data.

    Second dict may have extra keys, but must have all keys of first dict.
    """
    all_ok = True
    ok_exceptions = ['STAXOF', 'POLCALA', 'POLCALB', 'VELDEF', 'VELTYP','INSTRUME']

    for k in dict_a:
        if dict_b.has_key(k):
            try:
                if type(dict_a[k]) is float:
                    # Check that values are the same within tolerance for floats
                    assert np.allclose(dict_a[k], dict_b[k])
                elif type(dict_a[k]) is np.ndarray:
                    if type(dict_a[k][0]) is str:
                        assert all(dict_a[k] == dict_b[k])
                    else:
                        assert np.allclose(dict_a[k], dict_b[k])
                elif type(dict_a[k]) is np.core.defchararray.chararray:
                    assert all(dict_a[k] == dict_b[k])
                else:
                    assert dict_a[k] == dict_b[k]
                    assert type(dict_a[k]) == type(dict_b[k])
            except:
                if k not in ok_exceptions:
                    print "\nError:", k
                    print type(dict_a[k]), type(dict_b[k])
                    if type(dict_a[k]) is str and dict_a[k].strip() == '':
                        dict_a[k] = '(Empty str)'

                    if type(dict_b[k]) is str and dict_b[k].strip() == '':
                        dict_a[k] = '(Empty str)'

                    if type(dict_b[k]) in (float, int):
                         print dict_a[k], dict_b[k]
                    else:
                        print dict_a[k][:10], '\n', dict_b[k][:10]

                    try:
                        print "Len: %i %i"%(len(dict_a[k]), len(dict_b[k]))
                    except:
                        pass
                    all_ok = False
                else:
                    print "INFO: Known exception: %s"%k
        else:
            if k not in ok_exceptions:
                print "ERROR: %s not in both dictionaries"%k
                all_ok = False
            else:
                print "INFO: Known exception: %s"%k

    if all_ok:
        print "PASSED"
    else:
        print "ERROR"

    return all_ok