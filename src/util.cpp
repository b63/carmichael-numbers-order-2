#include <util.h>

/**
 * Convinience wrapper function for strchr.
 * Returns a pointer to the first occurence of `character`
 * in `str`. If not found, then pointer to end of null terminator
 * of null-terminated string `str` is returned.
 */
const char *
strchr_def(const char *str, int character)
{
    const char *rv { strchr(str, character) };
    if (rv) return rv;

    /* find end of str */
    for (rv=str; *rv; ++rv);
    return rv;
}


/**
 * Factorize `n` and append prime factors to `primes` and their correspoding powers to
 * `powers using trial division.` **Really slow**, useful only for one-time pre-computation.
 * Prime factors are of type NTL::ZZ so that even extremely large `n` can be factored.
 * Uses the NTL::NextPrime function.
 * @param  primes   vector to append primes factors to
 * @param  powers   vector to append correspoding power of the
 *                  prime factors to
 * @param  n        number  >= 1 to factorize. If n <= 1, then function
 *                  does nothing.
 */
void
factorize_slow(std::vector<NTL::ZZ> &primes, std::vector<long> &powers, const NTL::ZZ &n)
{
    if (n <= 1)
        return;

    Factorization f;
    NTL::ZZ p { 2 };
    NTL::ZZ m {n};
    NTL::ZZ sqrt {NTL::SqrRoot(m)};

    while (1)
    {
        long power {0};
        while (NTL::divide(m, m, p))
        {
            NTL::SqrRoot(sqrt, m);
            power++;
        }

        if(power)
        {
            primes.push_back(p);
            powers.push_back(power);
        }
        else if(p >= sqrt)
        {
            primes.push_back(m);
            powers.push_back(1);
            break;
        }

        if (NTL::IsOne(m)) break;
        NTL::NextPrime(p, p+1);
    }
}


std::unique_ptr<std::vector<size_t>>
random_subset(size_t set_size)
{
    static constexpr size_t bits {sizeof(unsigned long)<<3};
    static const unsigned long seed {(unsigned long)std::chrono::system_clock::now().time_since_epoch().count()};
    static std::mt19937_64 generator(seed);
    static std::uniform_int_distribution<unsigned long> dist(1, ULONG_MAX);

    std::unique_ptr<std::vector<size_t>> indicies { std::make_unique<std::vector<size_t>>() };
    while (indicies->size() == 0)
    {
        for(size_t k{0}; k < set_size; k += bits)
        {
            unsigned long rand_num {dist(generator)};
            while (rand_num)
            {
                int n {__builtin_ctzl(rand_num)};
                size_t index {n+k};
                if (index >= set_size)
                    break;

                indicies->push_back(index);
                rand_num &= rand_num-1;
            }
        }
    }

    return indicies;
}


/*
 * Returns the number of subsets of sizes ranging from 
 * `min_size` to `max_size` (both inclusive) when drawing
 * from a set of size `set_size`.
 */
size_t
calc_max_subsets(size_t set_size, size_t min_size, size_t max_size)
{
    size_t ret {0};
    max_size = max_size > 0 ? MIN(set_size, max_size) : set_size;
    for(; min_size <= max_size; min_size++)
        ret += binomial(set_size, min_size);
    return ret;
}


/*
 * Calculate the binomial coefficient C(n, m)
 * by constructing a pascal's triangle.
 */
size_t
binomial(unsigned long n, unsigned long m)
{
    if (m == 0 || n == m) return 1;
    else if (n == 0)      return 0;

    /* TODO: check the assembled output */
    /* allocate table  to store pascal's triangle */
    const size_t row {n+1};
    const size_t size {row*row};
    size_t *table {new size_t[size]};

    /* fill in pascal's triangle row by row */
    table[0] = 1;
    for(size_t i {1}; i <= n; ++i)
    {
        const size_t r {row*i};
        table[r] = i+1;
        table[r+i] = 1;
        /* TODO: only really need the last row, use
        *       some sort of cyclic data structure to
        *       reduce memory usage */
        for(size_t k{1}; k < i; ++k)
            table[r+k] = table[r-row+k] + table[r-row+k-1];
    }

    const size_t ret {table[row*(n-1)+m-1]};
    delete[] table;
    return ret;
}


/**
 * Parse a string of the form:
 *      "2^4 3^5  3^55 ... 10^5"
 *  as a `Factorization` object and return it
 */
Factorization
parse_factorization(const char *str)
{
    constexpr size_t BUF_SIZE { 50 };
    const char *prev { str };
    const char *ptr { strchr_def(str, ' ') };
    char buf[BUF_SIZE];

    long factor, power;
    Factorization f;

    while (prev && *prev)
    {
        if (*prev != ' ')
        {
            /* copy string to buffer */
            size_t diff { (size_t) MAX(0, ptr-prev) };
            strncpy(buf, prev, MIN(BUF_SIZE, diff+1));
            buf[MIN(BUF_SIZE-1,diff)] = '\0';

            parse_factor_string<long>(factor, power, buf);
            f.primes.push_back(factor);
            f.powers.push_back(power);
        }

        if (*ptr)
        {
            prev = ptr + 1;
            ptr = strchr_def(ptr+1, ' ');
        }
        else
        {
            prev = nullptr;
        }
    }
    /* TODO: check for NRVO */
    return f;
}


/**
 * Parse a string of space separated integers containing up 50 digits
 * and add them to `list`.
 * Returns the number of integers added to `list`.
 */
size_t
parse_numbers(std::vector<long> &list, const char *str)
{
    constexpr size_t BUF_SIZE { 51 };
    const char *prev { str };
    const char *ptr { strchr_def(str, ' ') };
    char buf[BUF_SIZE];

    size_t count { 0 };
    while (prev && *prev)
    {
        if (*prev != ' ')
        {
            /* copy string to buffer */
            size_t diff { (size_t) MAX(0, ptr-prev) };
            strncpy(buf, prev, MIN(BUF_SIZE, diff+1));
            buf[MIN(BUF_SIZE-1,diff)] = '\0';

	    list.push_back(atol(buf));
	    count++;
        }

        if (*ptr)
        {
            prev = ptr + 1;
            ptr = strchr_def(ptr+1, ' ');
        }
        else
        {
            prev = nullptr;
        }
    }

    return count;
}

int
parse_args(int argc, char **argv, long &max, Factorization &f)
{
    long factor, power;
    int limit { 0 };
    for (int i {1};  i < argc; i++)
    {
        if (*argv[i] == '-')
        {
            if (!limit && argv[i][1] == 'l')
            {
                size_t digits=0;
                char *ptr = argv[i]+2;
                for(; *ptr; digits++,ptr++)
                    if(*ptr < '0' || *ptr > '9') break;

                // at least 1 digits and no nondigit characters
                if (digits && !*ptr)
                {
                    max = NTL::conv<long>(argv[i]+2);
                    limit = 1;
                    continue;
                }
            }
        }
        parse_factor_string<long>(factor, power, argv[i]);
        f.primes.push_back(factor);
        f.powers.push_back(power);
    }

    if (f.primes.size() < 1)
        return 0;
    return 1;
}

/*
 * Ensures that `n` includes a `factor` as a factor.
 * Note: Assumes that the vectors in `n` are sorted from
 * least to greatest.
 * @param n       reference to Factorization object consisting of two sorted vectors
 *                that contains n's prime factors and their powers
 * @param factor  factorization of the factor that `n` must have
 */
Factorization
include_as_factor(const Factorization &n, const Factorization &factor)
{
    const size_t factor_size = factor.primes.size();
    const size_t n_size      = n.primes.size();

    std::unordered_map<long, long> map {};

    size_t count = 0;
    for(size_t i = 0; i < n_size; i++)
    {
        long power = n.powers[i];
        /* ignore primes with exponent 0 */
        if ( power && (map[n.primes[i]] += power) == power)
            count++; /* this prime hasn't been encountered before */
    }

    for(size_t i = 0; i < factor_size; i++)
    {
        long power = factor.powers[i];
        long &cur_power = map[factor.primes[i]];
        if (cur_power >= power)
            continue; /* ignore zero exponent primes & those with high enough exponent */

        if (!(cur_power = power)) 
            ++count;
    }

    Factorization new_n;
    new_n.primes.reserve(count);
    new_n.powers.reserve(count);

    // TODO: test if fater than building unordered list, then sorting using
    //       algorithm of nlog(n) complexity?
    for (auto it = map.cbegin(), cend = map.cend(); it != cend; ++it)
    {
        long prime = it->first, power = it->second;

        size_t i = 0, nn_size = new_n.primes.size();
        for (; i < nn_size; ++i)
            if (prime < new_n.primes[i])
                break;

        new_n.primes.insert(new_n.primes.cbegin()+i, prime);
        new_n.powers.insert(new_n.powers.cbegin()+i, power);
    }

    // TODO: check for copy-elision
    return new_n;
}

std::ostream&
operator<<(std::ostream &os, const Factorization &f)
{
    const size_t len = f.primes.size();
    if (len != f.powers.size())
    {
        throw std::length_error("primes and powers vectors have unequal lengths");
    }

    for (size_t k = 0; k < len; k++)
    {
        if (k > 0) os << " * ";
        os << f.primes[k] << "^" << f.powers[k];
    }

    return os;
}


// generate list of primes using a sieve and store
// them in `primes` vector.
// @param primes vector instance where the primes will be stored
// @param size   size of the sieve used to generate primes
//               effectively limits the largest primes.
//               The largest prime will be <= size
void sieve_primes(std::vector<long> &primes, const size_t size, 
        const std::function<bool(size_t, long)> *filter)
{
    // list of integers that will be sieved
    long *sieve_array = new long[size];
    for (size_t i = 0; i < size; i++) {
        sieve_array[i] = i+1;
    }

    sieve(sieve_array, size);

    // collect all the primes
    if (filter == nullptr)
    {
        for (size_t i = 1; i < size; i++)
        {
            long n = sieve_array[i];
            if (n != 0)
            {
                primes.push_back(n);
            }
        }
    }
    else
    {
        for (size_t i = 1; i < size; i++)
        {
            long n = sieve_array[i];
            if (n != 0 && (*filter)(i, n))
            {
                primes.push_back(n);
            }
        }
    }
}



// sieve an array of integers starting from 1 to create a  list of primes
// @prarm arr  pointer to the array, should contain [1, 2, ..., size-1]
// @param size the length of the array
void sieve(long *arr, size_t size)
{
    size_t index = 1;
    while (index < size)
    {
        long jump = arr[index];
        if (jump == 0)
        {
            index++;
            continue;
        }

        size_t  sieve_index = index + jump;

        for (; sieve_index < size; sieve_index += jump)
        {
            arr[sieve_index] = 0;
        }

        index++;
    }
}


bool poklington_test(const NTL::ZZ &N)
{
    static std::map<const NTL::ZZ, bool> _cache {};

    // for nonexistent key, _cache[key] is false
    // source: https://stackoverflow.com/questions/11058422/map-operator-and-bool-as-value
    if (_cache[N]) return true;

    const NTL::ZZ &N_1 {N-1};
    const NTL::ZZ &bound {NTL::SqrRoot(N)+1};

    // find f > sqrtN such that f | N-1
    NTL::ZZ f {1};
    NTL::ZZ r {N_1};
    std::vector<long> prime_factors {};
    for (size_t i {0}; f < bound; i++)
    {
        long ith_prime = get_nth_prime(i);
        // if ith_prime | r, then r /= ith_prime
        // otherwise r doesn't change
        if (NTL::divide(r, r, ith_prime))
        {
            prime_factors.push_back(ith_prime);
            f *= ith_prime;

            // keep dividing if possible
            while (f < bound && NTL::divide(r, r, ith_prime))
                f *= ith_prime;
        }
    }

    // TODO: need to ensure f and r are relatively prime
    return false;
}


/**
 * Returns the nth prime number using the PrimeSeq class
 * from the NTL library.
 */
long get_nth_prime(size_t n)
{
    static std::vector<long> _cache {0};
    static NTL::PrimeSeq s {};

    size_t size = _cache.size();
    if ( n < size)
    {
        return _cache[n];
    }

    _cache.reserve(n);
    while (size < n)
    {
        _cache.push_back(s.next());
        size++;
    }

    return _cache[n];
}


/** 
 * prints a carriage return ('\r') and `width` number
 * of spaces followed by another carriage return.
 * Helpful to 'erase' text written to a terminal
 */
void clrln(size_t width)
{
    char *spaces {new char[width+3]};
    spaces[0] = '\r';
    for (size_t i = 1; i <= width; i++) spaces[i] = ' ';
    spaces[width+1] = '\r';
    spaces[width+2] = '\0';
    std::cout << spaces;
    delete[] spaces;
}


/**
 * Multiplies each element of `factor` raised to the corresponding element of `powers`
 * and stores the final product in `prod`.
 *
 * @param prod    referece to NTL::ZZ where the final product is accumulated
 * @param factors vector of factors
 * @param powers  vector of powers for each factor in`factors`
 *
 * Returns a referce to `prod`.
 */
NTL::ZZ&
multiply_factors(NTL::ZZ &prod, const std::vector<long> &factors, const std::vector<long> &powers)
{
    prod = 1;
    const size_t size { factors.size() };
    for (size_t i = 0; i < size; i++)
    {
        NTL::ZZ pow { NTL::power_ZZ(factors[i], powers[i]) };
        prod *= pow;
    }

    return prod;
}


/**
 * returns the density as a float, i.e 2^|prime_size| / L
 */
NTL::RR
get_density(NTL::ZZ &L, size_t prime_size)
{
    static const NTL::RR log2 { NTL::log(NTL::conv<NTL::RR>(2)) };
    NTL::RR mag_P { NTL::conv<NTL::RR>(prime_size) };
    NTL::RR l { NTL::conv<NTL::RR>(L) };
    if ( l <= 0)
        return NTL::conv<NTL::RR>(0);
    else
        return mag_P - NTL::log(L)/log2;
}

