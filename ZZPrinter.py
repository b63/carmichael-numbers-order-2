import gdb
class ZZPrinter:
    "Basic pretty pinter for NTL:ZZ type"

    def __init__(self, val):
        self.val = val
        self.ptype = gdb.lookup_type('unsigned long').pointer()

    def to_string(self):
        _ntl_gbigint_ptr = self.val['rep']['rep']
        if _ntl_gbigint_ptr == 0:
            return '0'

        size = _ntl_gbigint_ptr.dereference()['size_']
        body =  (_ntl_gbigint_ptr+1).cast(self.ptype)

        sign = 1
        if size < 0:
            size = -size
            sign = -1

        longs = [int((body+i).dereference()) for i in range(size)]
        zz = 0
        base = 1
        for i in range(size):
            if i > 0:
                base *= 2**64
            zz += longs[i]*base

        return ('-' if sign < 0 else '') + str(zz)

def lookup_type(val):
    if str(val.type) == 'NTL::ZZ':
        return ZZPrinter(val)
    return None
gdb.pretty_printers.append(lookup_type)



