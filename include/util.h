#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <ostream>
#include <random>
#include <climits>
#include <chrono>
#include <stdio.h>
#include <cstdlib>
#include <climits>
#include <stdexcept>
#include <functional>
#include <map>
#include <unordered_map>
#include <memory>

#include <NTL/ZZ.h>
#include <NTL/RR.h>

#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))


/* Structures used for generating statistics at the end */
struct CoFactorSet
{
    /* number of primes associated with the set of cofactors */
    size_t num_primes;
    /* set of cofactors that have `num_primes` number of associated primes */
    std::vector<long> cofactors;
};


struct SizeCount
{
    /* number of primes */
    size_t num_primes;
    /* number of cofactors that have `num_primes` number of associated primes */
    size_t num_cofactors;
};

struct Factorization
{
    std::vector<long> primes;
    std::vector<long> powers;
};

std::ostream& operator<<(std::ostream& os, const Factorization& f);

/*************** FUNCTIONS ***********/

// clear `width` number of characters from terminal
// by writing that number of spaces followed by a
// carriage return
void clrln(size_t width=80);

// stores list of primes < size in `primes` vector
void sieve_primes(std::vector<long> &primes, const size_t size, 
        const std::function<bool(size_t, long)> *filter = nullptr);
// sieves `arr` , and `size` is lenght of arr
void sieve(long *arr, size_t size);

// uses NTL's built in prime generator to
// generate consecutive primes. Counts start from 0,
// so get_nth_prime(0) returns 2.
long get_nth_prime(size_t n);

NTL::ZZ& multiply_factors(NTL::ZZ &prod, const std::vector<long> &factors, const std::vector<long> &powers);

NTL::RR get_density(NTL::ZZ &L, size_t prime_size);

Factorization include_as_factor(const Factorization &n, const Factorization &factor);

const char * strchr_def(const char *str, int character);

Factorization parse_factorization(const char *str);

size_t parse_numbers(std::vector<long> &list, const char *str);

int parse_args(int argc, char **argv, long &max, Factorization &f);

size_t calc_max_subsets(size_t set_size, size_t min_size, size_t max_size);

size_t binomial(unsigned long n, unsigned long m);

std::unique_ptr<std::vector<size_t>> random_subset(size_t set_size);

void factorize_slow(std::vector<NTL::ZZ> &primes, std::vector<long> &powers, const NTL::ZZ &n);

NTL::ZZ eulers_toitent(const NTL::ZZ &n);

size_t estimate_subsets_size_bool(size_t set_size, size_t min_size, size_t max_size);

std::vector<std::vector<long>>
split_vector(const std::vector<long> &src, size_t n, bool shuffle=false);

/**************** TEMPLATE FUNCTIONS **************************/

template <typename T>
std::vector<T>
join_partitions(const std::vector<std::vector<T>> &partitions, size_t start, size_t end)
{
    std::vector<T> joined;
    size_t total_size {0};
    for (size_t i = start; i < end; i++)
        total_size += partitions[i].size();

    joined.reserve(total_size);
    for (size_t i = start; i < end; i++)
    {
        for (const auto &elem : partitions[i])
        {
            joined.push_back(elem);
        }
    }

    return joined;
}



template <typename T>
std::vector<T>
join_vector(const std::vector<T> &v0, const std::vector<T> &v1)
{
    std::vector<T> joined;
    joined.reserve(v0.size() + v1.size());
    for(const auto &x : v0) joined.push_back(x);
    for(const auto &x : v1) joined.push_back(x);

    return joined;
}


template <typename T>
inline T bound(const T &val, const T &min, const T &max)
{
    if (val < min)
        return min;
    else if (val > max)
        return max;
    else
        return val;
}


template <typename T, typename V>
V product(T iterable)
{
    V ret {1};
    for (auto v : iterable)
        ret *= v;
    return ret;
}


/**
 * Get the lcm of a container of ZZ's.
 * Must have at least one element.
 */
template <typename T>
NTL::ZZ get_lcm(const T &arr, size_t size)
{
    if (size < 2)
        return arr[0];
    NTL::ZZ gcd {NTL::GCD(arr[0], arr[1])};
    NTL::ZZ prod {arr[0] * arr[1]};
    NTL::ZZ lcm {prod/gcd};

    for(size_t i=2; i < size; ++i)
    {
        gcd = NTL::GCD(lcm, arr[i]);
        prod = lcm * arr[i];
        lcm = prod / gcd;
    }

    return lcm;
}


/**
 * Splits the elements of vector<T> `set` into two halves.
 * If n is the number of elements in `set`, the first
 * floor(n/2) elements will be appended to `first`, and
 * the rest will be appended to `second`.
 */
template <typename T>
void split_half(std::vector<T> &first, std::vector<T> &second,
        const std::vector<T> set)
{
    const size_t set_size {set.size()};
    size_t first_size {first.size()};
    first.reserve(first_size + set_size/2);
    for(auto it{set.cbegin()},end=it+set_size/2; it != end; ++it)
        first.push_back(*it);

    first_size = first.size() - first_size;
    second.reserve(second.size() + set_size-first_size);
    for(auto it{set.cbegin()+first_size},end=set.cend(); it != end; ++it)
        second.push_back(*it);
}



/*
 * Returns the number of bytes storing subsets of sizes ranging from
 * `min_size` to `max_size` (both inclusive) would take when drawing
 * from a set of size `set_size` if each element takes sizeof(T) bytes,
 * where T is the template parameter.
 */
template <typename T>
size_t estimate_subsets_size(size_t set_size, size_t min_size, size_t max_size)
{
    const size_t per_T {sizeof(T)};
    size_t ret {0};
    max_size = max_size > 0 ? MIN(set_size, max_size) : set_size;
    for(; min_size <= max_size; min_size++)
        ret += per_T * binomial(set_size, min_size) * min_size;
    return ret;
}



// prints the integers in arr specified by the array of indices ind 
// with a '*' as separator
template <typename T>
void printProd(const std::vector<size_t> &ind, const std::vector<T> &arr)
{
    size_t factors = ind.size();
    for (size_t k = 0; k < factors; k++)
    {
        if (k > 0) std::cout << " * ";
        std::cout << arr[ind[k]];
    }
}


/**
 * prints a vector of type T as a formatted list
 *      {arr[0], arr[1], .., arr[n-1]}
 */
template <typename T>
void printVec(const std::vector<T> &arr) 
{
    const size_t size = arr.size();
    std::cout << "{";

    for (size_t i = 0; i < size; i++)
    {
        if (i > 0) 
        {
            std::cout << ", ";
        }
        std::cout << arr[i];
    }

    std::cout << "}";
}


template <typename T>
void print_bool_vec(const std::vector<bool> &arr, const std::vector<T> &src)
{
    const size_t size = arr.size() < src.size() ? arr.size() : src.size();
    std::cout << "{";
    bool first = true;

    for (size_t i = 0; i < size; i++)
    {
        if (arr[i])
        {
            if (first) first = false;
            else       std::cout << ", ";
            std::cout << src[i];
        }
    }

    std::cout << "}";
}


template <typename T>
std::ostream& operator<<(std::ostream &os, const std::vector<T> &arr) 
{
    size_t size = arr.size();
    os << "{";

    for (size_t i = 0; i < size; i++)
    {
        if (i > 0) 
        {
            os << ", ";
        }
        os << arr[i];
    }

    os << "}";
    return os;
}


template <typename T, long unsigned int N>
std::ostream& operator<<(std::ostream &os, const std::array<T,N> &arr) 
{
    size_t size = arr.size();
    os << "{";

    for (size_t i = 0; i < size; i++)
    {
        if (i > 0) 
        {
            os << ", ";
        }
        os << arr[i];
    }

    os << "}";
    return os;
}


/**
 * Prints the prime factorization i.e "2^6 * 3^2  * 7^5"
 * @param factors vector of factors
 * @param powers  vector of corresponding powers 
 */
template <typename T>
void printFactorization(const std::vector<T> &factors, const std::vector<T> &powers)
{
    size_t len = factors .size();
    for (size_t k = 0; k < len; k++)
    {
        if (k > 0) std::cout << " * ";
        std::cout << factors[k] << "^" << powers[k];
    }
}


/**
 * Binary searches elements of a collection using randomly acessible iterator
 * starting at offset `start` and the next `length` elements.
 * The element that is being searched for is the element for which `comp`
 * function returns 0.
 *
 * If the element `x` being searched is compared to some other element `y`, then
 * `comp` should return;
 *      0    if  x == y
 *      <0    if  x < y
 *      >0    if  x > y
 * The offset where the element being serached for is found is written to `ind`.
 * If the element is not found, the the offset where the element would be found
 * if it exsisted is written to `ind`. Meaning, if the element being searched
 * for is inserted to the collected at the offset written to `ind`, the order
 * of the collection is perserved.
 *
 * Returns true if element was found, otherwise false.
 *
 * Ex. to binary search a vector<int> vec for `7`:
 *       binary_search(&ind, vec.begin(), 0, vec.size(), [](int y){ return x-y; } )
 *
 * @param ind     reference to variable where the offset is written
 * @param it      random access iterator
 * @param start   offset to start the search
 * @param length  number of elments from `start` to search
 * @param comp    function that takes an arugment and returns and integer whose
 *                sign determines whether the element being searched for 
 *                is less than or greater than the argument
 */
template <typename RandomAccessIterator, typename T>
bool
binary_search(
        size_t &ind, const RandomAccessIterator it, size_t start, size_t length,
        std::function<int(const T&)> comp)
{
    if (start == length) {
        ind = start;
        return false;
    }

    size_t mid;
    int order;
    do 
    {
        mid = (start + length)/2;
        order = comp(it[mid]);

        if (order < 0)
        {
            start = mid+1;
        }
        else if (order > 0)
        {
            length = mid;
        }
        else
        {
            ind = mid;
            return true;
        }
    } while (start < length);

    ind = mid < start ? start : mid;
    return false;
}


template <typename T>
void 
print_stats(const std::map<const long, std::vector<T>> &factor_map, const NTL::ZZ &L_prime, 
        size_t max_ranking, size_t num_print_k = 5)
{

    std::vector<CoFactorSet> ranking;
    std::vector<SizeCount> bag_sizes;

    for (auto it=factor_map.cbegin(); it != factor_map.cend(); it++)
    {
        const long k      { it->first         };
        const size_t size { it->second.size() };

        /* update bag_sizes */
        size_t bag_size { bag_sizes.size() };
        size_t ind {0};
        bool found { 
            binary_search<std::vector<SizeCount>::const_iterator, SizeCount>
                (
                    ind, bag_sizes.begin(), 0, bag_size, 
                    [=](const SizeCount &item)->int 
                    {
                        if      (item.num_primes < size) return -1;
                        else if (item.num_primes > size) return  1;
                        else                             return  0;
                    }
                )
            };

        if (found)
        {
            bag_sizes[ind].num_cofactors += 1;
        }
        else
        {
            SizeCount entry {size, 1};
            if (ind == bag_size)
                bag_sizes.push_back(entry);
            else
                bag_sizes.insert(bag_sizes.begin()+ind, entry);
        }

        /* update ranking */
        found = binary_search<std::vector<CoFactorSet>::const_iterator, CoFactorSet>
            (
                ind, ranking.begin(), 0, ranking.size(),
                [=](const CoFactorSet &item)->int
                {
                    if      (item.num_primes < size) return -1;
                    else if (item.num_primes > size) return  1;
                    else                             return  0;
                }
            );
        if (found)
        {
            ranking[ind].cofactors.push_back(k);
        }
        else
        {
            CoFactorSet entry {size, std::vector<long>{k}};
            ranking.insert(ranking.begin()+ind, entry);
            if (ranking.size() > max_ranking)
                ranking.erase(ranking.begin());
        }

    }

    if (ranking.size() > 0)
    {
        const CoFactorSet &c {ranking[ranking.size()-1]};
        if (c.cofactors.size() > 0)
        {
            NTL::ZZ k {c.cofactors[0]};
            NTL::ZZ L {L_prime * k};
            NTL::RR density;
            try
            {
                density = get_density(L, c.num_primes);
                std::cout << "density = " << density << "\n";
                std::cout << "L = " << L << "\n";
                std::cout << "  = " << L_prime << " * " << k << "\n";
            }
            catch (std::exception &e)
            {
                std::cout << "Exception: " << e.what() << "\n";
            }
        }
    }

    /* print bag */
    size_t w {10};
    std::cout << std::setw(w) << "# primes" << "," << std::setw(w) << "# co-factors" << "\n";
    for(size_t i = 0; i < bag_sizes.size(); i++)
    {
        const SizeCount &item {bag_sizes[bag_sizes.size() - i - 1]};
        std::cout << std::setw(w) << item.num_primes << ",";
        std::cout << std::setw(w) << item.num_cofactors << "\n";
    }
    std::cout << "\n";

    /* print ranking */
    if (ranking.size() > 0)
    {
        w = 5;
        for (auto it {ranking.crbegin()}; it != ranking.crend(); it++)
        {
            std::cout << "num primes=" << it->num_primes << "\n";
            const std::vector<long> &cfs {it->cofactors};
            size_t j = 0;
            for (; j < num_print_k && j < cfs.size(); j++)
            {
                std::cout << "k=" << std::setw(w) << cfs[j] << " -> ";
                const std::vector<T> &primes {factor_map.at(cfs[j])};
                printVec<T>(primes);
                std::cout << "\n";
            }

            if (j < cfs.size())
                std::cout << "(" << (cfs.size()-j) << " more...)\n";

            if (it+1 < ranking.crend())
                std::cout << "\n";
        }
    }
}


/**
 * Parses null-terminated `str` as "factor^power", and writes them to 
 * `factor` and `power`. If 'power' is not given, then 1 is assumed.
 * Throws std::invalid_argument if there were parsing errors.
 * Returns true if factor and power were successfully extracted.
 *
 * @param factor  reference to variable where factor is to be written
 * @param power   reference to variable where the power is to be written
 * @param str     pointer to first character of string to parse
 */
template <typename T>
bool
parse_factor_string(T &factor, T &power, char *str)
{
    size_t i { 0 };
    char *ch {str}, *mid {0};
    while(*ch)
    {
        if (*ch < '0' || *ch > '9')
        {
            /* character is not a digit */
            if (!i)
            {
                char buf[100];
                snprintf(buf, 100, "first character cannot be non-digit, '%s'", str);
                throw std::invalid_argument(buf);
            }
            else if (mid)
            {
                char buf[100];
                snprintf(buf, 100, "only digit characters allowed, '%s'", str);
                throw std::invalid_argument(buf);
            }
            else
            {
                mid = ch+1;
            }
        }
        i++;
        ch++;
    }
    if(!i) return false;

    if (!mid)
    {
        power = NTL::conv<T>(1);
        factor = NTL::conv<T>(str);
    }
    else
    {
        power = NTL::conv<T>(mid);
        char prev { *mid };
        *mid = 0;
        *mid = prev;
        factor = NTL::conv<T>(str);
    }

    return true;
}


#endif
