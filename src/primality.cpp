#include <vector>
#include <NTL/ZZ.h>

#include <primality.h>

/* ZZ version of Pocklington test.  This algorithm attempts to prove a number is prime.
If the number is prime, the running time is expected polynomial.  If the number is composite,
it may be very slow, so make sure to apply a primality test first.
Input is n along with a factorization of n-1 given by a pair of vector<longs>
Output is a bool.  True means proven prime, False means failed to prove prime.
*/
bool
pocklington(const NTL::ZZ &n, const std::vector<long> &primes)
{
    /* PowerMod will raise error if given 1 as an argument */
    if (n == 2) return true;

    // in English: for every prime q | n-1, find a with a^(n-1) = 1 mod n, gcd(a^{(n-1)/q}-1, n) =1 
    // first get (pointers to) the primes of f

    // Output boolean is True unless proven otherwise
    bool output = true;

    // boolean for the inner loop
    bool done;
    // these are used in the inner loop for powers needed
    NTL::ZZ full_power;
    NTL::ZZ sans_q_power;
    // these variables will store gcds
    NTL::ZZ g1; NTL::ZZ g2; NTL::ZZ g3;

    // Outer loop is over the prime_factors vector, so checking all q | n-1
    for (size_t i {0}; i < primes.size(); i++) {
        //cout << "prime = " << ps->at(i) << "\n";

        // break out of the loop if output is false, since that means n proven composite
        if (output == false) {
            break;
        }

        // we are looking for a generator a.  We start at 2 and continue until done
        done = false;
        NTL::ZZ a {2};
        while (!done) {
            //cout << "a = " << a << "\n";

            // compute the two powers needed
            NTL::PowerMod(full_power, a, n-1, n);
            NTL::PowerMod(sans_q_power, a, (n-1)/primes[i], n);
            sans_q_power -= 1;

            // compute gcds with n.  Unlikely to catch a factor, but worth checking
            NTL::GCD(g1, a, n);
            NTL::GCD(g2, full_power, n);
            NTL::GCD(g3, sans_q_power, n);

            // if gcd(a,n) is nontrivial, then n composite
            if (!NTL::IsOne(g1) && g1 != n) {
                output = false;
                done = true;
            }
            // if gcd(full_power, n) is nontrivial, then n composite
            else if (!NTL::IsOne(g2) && g2 != n) {
                output = false;
                done = true;
            }
            // if gcd(sans_q_power, n) is nontrivial, then n composite
            else if (!NTL::IsOne(g3) && g3 != n) {
                output = false;
                done = true;
            }
            // if a^(n-1) = 1 mod n and gcd(a^{(n-1)/q}-1, n) = 1, we've passed for this prime
            else if (NTL::IsOne(full_power) && NTL::IsOne(g3)) {
                done = true;
                //cout << "for q = " << ps->at(i) << " found a = " << a << "\n";
            }
            // otherwise increment a and try again.  Stop if a = n
            else {
                a++;
                if (a == n) {
                    done = true;
                    output = false;
                }
            }
        } // end while loop

    } // end for loop
    return output;
}
