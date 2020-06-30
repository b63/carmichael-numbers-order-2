weights = [1, 1.2, 2.2, 3.4]

def get_inc_possible(weight_i, value_distribution_set, cur_range=None, bound=None):
    """
    weight_i                       index into weights list
    cur_range                      tuple(low, high)
    value_distribution_set         list of tuples (value, dict)
    """
    if cur_range is None:
        cur_range = (0, weights[weight_i])

    low, high = cur_range

    new_set = []
    weight = weights[weight_i]
    nlow, nhigh = low + weight, high + weight
    if bound is not None:
        nlow  = min(nlow, bound[0])
        nhigh = min(nhigh, bound[1])


    for value, distribution in value_distribution_set:
        i = 0;
        while value < nlow:
            value += weight
            i += 1

        weight_inc = weight_i in distribution
        while nlow <= value < nhigh:
            copy = distribution.copy()
            if weight_inc:
                copy[weight_i] += i
            elif i > 0:
                copy[weight_i] = i
            new_set.append( (value, copy) )
            value += weight
            i += 1

    return (nlow, nhigh), new_set


def get_all(bound):
    low, high = bound

    cur_set = [(0, {})]
    ranges_set = []

    for i, weight in enumerate(weights):
        new_range, new_subset = get_inc_possible(i, cur_set, cur_range=None, bound=(low,high))
        if len(new_subset) > 0:
            ranges_set.append((new_range, new_subset))

    nweights_i = [i for i in range(len(weights))]
    while len(nweights_i) > 0:
        nset = []
        j = 0
        while j < len(nweights_i): 
            i = nweights_i[j]
            quit = True
            for r, s in ranges_set:
                nr, ns = get_inc_possible(i, s, cur_range=r, bound=(low,high))
                #print("{} {} -> {}, {}".format(str(r), str(s), str(nr), str(ns)))
                nset.append((nr, ns))
                if nr != (low,high):
                    quit = False

            if quit:
                nweights_i.pop(j)
            else:
                j += 1
        ranges_set = nset

    dset = set()
    for r, s in ranges_set:
        for v, d in s:
            fs = frozenset(d.items())
            dset.add(fs)

    print(weights)
    for d in dset:
        arr = [0 for i in range(len(weights))]
        for i, v in d:
            arr[i] = v
        print(arr)

    return dset











