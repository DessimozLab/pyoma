from future.builtins import str
import numpy
import re
import tables

"""This package is intended to parse the darwin locus structure
and create a numpy recarray out of it. It is put in an own
package to protect the global namespace from being polluted with
the global variable."""

eNr = -1

def parse(omaEntryNr, s):
    global eNr
    transInput = re.sub(r'(Before\(|After\()?(?P<r1>\d+)(?(1)\))\s*\.\.\s*'
                         '(:Before\(|After\()?(?P<r2>\d+)(?(3)\))',
                        r'(\g<r1>,\g<r2>)', s)
    eNr = omaEntryNr
    try:
        res = eval(transInput)
        if not isinstance(res, numpy.recarray):
            res = join(res)
        return res
    except Exception as e:
        raise ValueError('cannot parse "'+transInput+'" for entry '+str(eNr))


def join(*args):
    from .convert import LocusTable
    tab = numpy.rec.array(None, shape=(len(args),), 
        dtype=tables.file.dtype_from_descr(LocusTable))
    for i, op in enumerate(args):
        if len(op)==2:
            tab[i] =(eNr, op[0], op[1], 1)
        elif len(op)==3:
            tab[i] = (eNr, op[0], op[1], op[2])
        else:
            raise ValueError('unexpected length for '+op)
    return tab


def complement(op):
    if isinstance(op, numpy.recarray):
        op['Strand'] *= -1
        return op
    elif len(op) != 2:
        raise ValueError('unexpeted input for complement')
    return (op[0],op[1],-1,)

def Before(op):
    return op

def After(op):
    return op

