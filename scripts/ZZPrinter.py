import gdb
class ZZPrinter:
    "Basic pretty pinter for NTL:ZZ and bigint type"

    def __init__(self, val):
        self.val = val
        self.ptype = gdb.lookup_type('unsigned long').pointer()

    def parse_bigint(self, _ntl_gbigint_ptr):
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

    def to_string(self):
        t = str(self.val.type)
        if 'NTL::ZZ' in t:
            _ntl_gbigint_ptr = self.val['rep']['rep']
        elif '_ntl_gbigint' in t:
            _ntl_gbigint_ptr = self.val
        else:
            _ntl_gbigint_ptr = 0

        if _ntl_gbigint_ptr == 0:
            return '0'

        return self.parse_bigint(_ntl_gbigint_ptr)


def lookup_type(val):
    t = str(val.type)
    if 'NTL::ZZ' in t or '_ntl_gbigint' in t:
        return ZZPrinter(val)
    return None
gdb.pretty_printers.append(lookup_type)



