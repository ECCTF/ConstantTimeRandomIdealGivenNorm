#include <quaternion.h>
#include "internal.h"

/** @file
 *
 * @authors Antonin Leroux
 *
 * @brief Functions related to norm equation solving or special extremal orders
 */

void
quat_lattice_O0_set(quat_lattice_t *O0)
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_set(&(O0->basis[i][j]), 0);
        }
    }
    ibz_set(&(O0->denom), 2);
    ibz_set(&(O0->basis[0][0]), 2);
    ibz_set(&(O0->basis[1][1]), 2);
    ibz_set(&(O0->basis[2][2]), 1);
    ibz_set(&(O0->basis[1][2]), 1);
    ibz_set(&(O0->basis[3][3]), 1);
    ibz_set(&(O0->basis[0][3]), 1);
}

void
quat_lattice_O0_set_extremal(quat_p_extremal_maximal_order_t *O0)
{
    ibz_set(&O0->z.coord[1], 1);
    ibz_set(&O0->t.coord[2], 1);
    ibz_set(&O0->z.denom, 1);
    ibz_set(&O0->t.denom, 1);
    O0->q = 1;
    quat_lattice_O0_set(&(O0->order));
}

void
quat_order_elem_create(quat_alg_elem_t *elem,
                       const quat_p_extremal_maximal_order_t *order,
                       const ibz_vec_4_t *coeffs,
                       const quat_alg_t *Bpoo)
{

    // var dec
    quat_alg_elem_t quat_temp;

    // var init
    quat_alg_elem_init(&quat_temp);

    // elem = x
    quat_alg_scalar(elem, &(*coeffs)[0], &ibz_const_one);

    // quat_temp = i*y
    quat_alg_scalar(&quat_temp, &((*coeffs)[1]), &ibz_const_one);
    quat_alg_mul(&quat_temp, &order->z, &quat_temp, Bpoo);

    // elem = x + i*y
    quat_alg_add(elem, elem, &quat_temp);

    // quat_temp = z * j
    quat_alg_scalar(&quat_temp, &(*coeffs)[2], &ibz_const_one);
    quat_alg_mul(&quat_temp, &order->t, &quat_temp, Bpoo);

    // elem = x + i* + z*j
    quat_alg_add(elem, elem, &quat_temp);

    // quat_temp = t * j * i
    quat_alg_scalar(&quat_temp, &(*coeffs)[3], &ibz_const_one);
    quat_alg_mul(&quat_temp, &order->t, &quat_temp, Bpoo);
    quat_alg_mul(&quat_temp, &quat_temp, &order->z, Bpoo);

    // elem =  x + i*y + j*z + j*i*t
    quat_alg_add(elem, elem, &quat_temp);

    quat_alg_elem_finalize(&quat_temp);
}

int
quat_represent_integer(quat_alg_elem_t *gamma,
                       const ibz_t *n_gamma,
                       int non_diag,
                       const quat_represent_integer_params_t *params)
{

    if (ibz_is_even(n_gamma)) {
        return 0;
    }
    // var dec
    int found;
    ibz_t cornacchia_target;
    ibz_t adjusted_n_gamma, q;
    ibz_t bound, sq_bound, temp;
    ibz_t test;
    ibz_vec_4_t coeffs; // coeffs = [x,y,z,t]
    quat_alg_elem_t quat_temp;

    if (non_diag)
        assert(params->order->q % 4 == 1);

    // var init
    found = 0;
    ibz_init(&bound);
    ibz_init(&test);
    ibz_init(&temp);
    ibz_init(&q);
    ibz_init(&sq_bound);
    ibz_vec_4_init(&coeffs);
    quat_alg_elem_init(&quat_temp);
    ibz_init(&adjusted_n_gamma);
    ibz_init(&cornacchia_target);

    ibz_set(&q, params->order->q);

    // this could be removed in the current state
    int standard_order = (params->order->q == 1);

    // adjusting the norm of gamma (multiplying by 4 to find a solution in an order of odd level)
    if (non_diag || standard_order) {
        ibz_mul(&adjusted_n_gamma, n_gamma, &ibz_const_two);
        ibz_mul(&adjusted_n_gamma, &adjusted_n_gamma, &ibz_const_two);
    } else {
        ibz_copy(&adjusted_n_gamma, n_gamma);
    }
    // computation of the first bound = sqrt (adjust_n_gamma / p - q)
    ibz_div(&sq_bound, &bound, &adjusted_n_gamma, &((params->algebra)->p));
    ibz_set(&temp, params->order->q);
    ibz_sub(&sq_bound, &sq_bound, &temp);
    ibz_sqrt_floor(&bound, &sq_bound);

    // the size of the search space is roughly n_gamma / (p√q)
    ibz_t counter;
    ibz_init(&counter);
    ibz_mul(&temp, &temp, &((params->algebra)->p));
    ibz_mul(&temp, &temp, &((params->algebra)->p));
    ibz_sqrt_floor(&temp, &temp);
    ibz_div(&counter, &temp, &adjusted_n_gamma, &temp);

    // entering the main loop
    while (!found && ibz_cmp(&counter, &ibz_const_zero) != 0) {
        // decreasing the counter
        ibz_sub(&counter, &counter, &ibz_const_one);

        // we start by sampling the first coordinate
        ibz_rand_interval(&coeffs[2], &ibz_const_one, &bound);

        // then, we sample the second coordinate
        // computing the second bound in temp as sqrt( (adjust_n_gamma - p*coeffs[2]²)/qp )
        ibz_mul(&cornacchia_target, &coeffs[2], &coeffs[2]);
        ibz_mul(&temp, &cornacchia_target, &(params->algebra->p));
        ibz_sub(&temp, &adjusted_n_gamma, &temp);
        ibz_mul(&sq_bound, &q, &(params->algebra->p));
        ibz_div(&temp, &sq_bound, &temp, &sq_bound);
        ibz_sqrt_floor(&temp, &temp);

        if (ibz_cmp(&temp, &ibz_const_zero) == 0) {
            continue;
        }
        // sampling the second value
        ibz_rand_interval(&coeffs[3], &ibz_const_one, &temp);

        // compute cornacchia_target = n_gamma - p * (z² + q*t²)
        ibz_mul(&temp, &coeffs[3], &coeffs[3]);
        ibz_mul(&temp, &q, &temp);
        ibz_add(&cornacchia_target, &cornacchia_target, &temp);
        ibz_mul(&cornacchia_target, &cornacchia_target, &((params->algebra)->p));
        ibz_sub(&cornacchia_target, &adjusted_n_gamma, &cornacchia_target);
        assert(ibz_cmp(&cornacchia_target, &ibz_const_zero) > 0);

        // applying cornacchia
        if (ibz_probab_prime(&cornacchia_target, params->primality_test_iterations))
            found = ibz_cornacchia_prime(&(coeffs[0]), &(coeffs[1]), &q, &cornacchia_target);
        else
            found = 0;

        if (found && non_diag && standard_order) {
            // check that we can divide by two at least once
            // the treatmeat depends if the basis contains (1+j)/2 or (1+k)/2
            // we must have x = t mod 2 and y = z mod 2
            // if q=1 we can simply swap x and y
            if (ibz_is_odd(&coeffs[0]) != ibz_is_odd(&coeffs[3])) {
                ibz_swap(&coeffs[1], &coeffs[0]);
            }
            // we further check that (x-t)/2 = 1 mod 2 and (y-z)/2 = 1 mod 2 to ensure that the
            // resulting endomorphism will behave well for dim 2 computations
            found = found && ((ibz_get(&coeffs[0]) - ibz_get(&coeffs[3])) % 4 == 2) &&
                    ((ibz_get(&coeffs[1]) - ibz_get(&coeffs[2])) % 4 == 2);
        }
        if (found) {

#ifndef NDEBUG
            ibz_set(&temp, (params->order->q));
            ibz_mul(&temp, &temp, &(coeffs[1]));
            ibz_mul(&temp, &temp, &(coeffs[1]));
            ibz_mul(&test, &(coeffs[0]), &(coeffs[0]));
            ibz_add(&temp, &temp, &test);
            assert(0 == ibz_cmp(&temp, &cornacchia_target));

            ibz_mul(&cornacchia_target, &(coeffs[3]), &(coeffs[3]));
            ibz_mul(&cornacchia_target, &cornacchia_target, &(params->algebra->p));
            ibz_mul(&temp, &(coeffs[1]), &(coeffs[1]));
            ibz_add(&cornacchia_target, &cornacchia_target, &temp);
            ibz_set(&temp, (params->order->q));
            ibz_mul(&cornacchia_target, &cornacchia_target, &temp);
            ibz_mul(&temp, &(coeffs[0]), &coeffs[0]);
            ibz_add(&cornacchia_target, &cornacchia_target, &temp);
            ibz_mul(&temp, &(coeffs[2]), &coeffs[2]);
            ibz_mul(&temp, &temp, &(params->algebra->p));
            ibz_add(&cornacchia_target, &cornacchia_target, &temp);
            assert(0 == ibz_cmp(&cornacchia_target, &adjusted_n_gamma));
#endif
            // translate x,y,z,t into the quaternion element gamma
            quat_order_elem_create(gamma, (params->order), &coeffs, (params->algebra));
#ifndef NDEBUG
            quat_alg_norm(&temp, &(coeffs[0]), gamma, (params->algebra));
            assert(ibz_is_one(&(coeffs[0])));
            assert(0 == ibz_cmp(&temp, &adjusted_n_gamma));
            assert(quat_lattice_contains(NULL, &((params->order)->order), gamma));
#endif
            // making gamma primitive
            // coeffs contains the coefficients of primitivized gamma in the basis of order
            quat_alg_make_primitive(&coeffs, &temp, gamma, &((params->order)->order));

            if (non_diag || standard_order)
                found = (ibz_cmp(&temp, &ibz_const_two) == 0);
            else
                found = (ibz_cmp(&temp, &ibz_const_one) == 0);
        }
    }

    if (found) {
        // new gamma
        ibz_mat_4x4_eval(&coeffs, &(((params->order)->order).basis), &coeffs);
        ibz_copy(&gamma->coord[0], &coeffs[0]);
        ibz_copy(&gamma->coord[1], &coeffs[1]);
        ibz_copy(&gamma->coord[2], &coeffs[2]);
        ibz_copy(&gamma->coord[3], &coeffs[3]);
        ibz_copy(&gamma->denom, &(((params->order)->order).denom));
    }
    // var finalize
    ibz_finalize(&counter);
    ibz_finalize(&bound);
    ibz_finalize(&temp);
    ibz_finalize(&sq_bound);
    ibz_vec_4_finalize(&coeffs);
    quat_alg_elem_finalize(&quat_temp);
    ibz_finalize(&adjusted_n_gamma);
    ibz_finalize(&cornacchia_target);
    ibz_finalize(&q);
    ibz_finalize(&test);

    return found;
}

int
quat_sampling_random_ideal_O0_given_norm(quat_left_ideal_t *lideal,
                                         const ibz_t *norm,
                                         int is_prime,
                                         const quat_represent_integer_params_t *params,
                                         const ibz_t *prime_cofactor)
{

    ibz_t n_temp, norm_d;
    ibz_t disc;
    quat_alg_elem_t gen, gen_rerand;
    int found = 0;
    ibz_init(&n_temp);
    ibz_init(&norm_d);
    ibz_init(&disc);
    quat_alg_elem_init(&gen);
    quat_alg_elem_init(&gen_rerand);

    // when the norm is prime we can be quite efficient
    // by avoiding to run represent integer
    // the first step is to generate one ideal of the correct norm
    if (is_prime) {

        // we find a quaternion element of norm divisible by norm
        while (!found) {
            // generating a trace-zero element at random
            ibz_set(&gen.coord[0], 0);
            ibz_sub(&n_temp, norm, &ibz_const_one);
            for (int i = 1; i < 4; i++)
                ibz_rand_interval(&gen.coord[i], &ibz_const_zero, &n_temp);

            // first, we compute the norm of the gen
            quat_alg_norm(&n_temp, &norm_d, &gen, (params->algebra));
            assert(ibz_is_one(&norm_d));

            // and finally the negation mod norm
            ibz_neg(&disc, &n_temp);
            ibz_mod(&disc, &disc, norm);
            // now we check that -n is a square mod norm
            // and if the square root exists we compute it
            found = ibz_sqrt_mod_p(&gen.coord[0], &disc, norm);
            found = found && !quat_alg_elem_is_zero(&gen);
        }
    } else {
        assert(prime_cofactor != NULL);
        // if it is not prime or we don't know if it is prime, we may just use represent integer
        // and use a precomputed prime as cofactor
        assert(!ibz_is_zero(norm));
        ibz_mul(&n_temp, prime_cofactor, norm);
        found = quat_represent_integer(&gen, &n_temp, 0, params);
        found = found && !quat_alg_elem_is_zero(&gen);
    }
#ifndef NDEBUG
    if (found) {
        // first, we compute the norm of the gen
        quat_alg_norm(&n_temp, &norm_d, &gen, (params->algebra));
        assert(ibz_is_one(&norm_d));
        ibz_mod(&n_temp, &n_temp, norm);
        assert(ibz_cmp(&n_temp, &ibz_const_zero) == 0);
    }
#endif

    // now we just have to rerandomize the class of the ideal generated by gen
    found = 0;
    while (!found) {
        for (int i = 0; i < 4; i++) {
            ibz_rand_interval(&gen_rerand.coord[i], &ibz_const_one, norm);
        }
        quat_alg_norm(&n_temp, &norm_d, &gen_rerand, (params->algebra));
        assert(ibz_is_one(&norm_d));
        ibz_gcd(&disc, &n_temp, norm);
        found = ibz_is_one(&disc);
        found = found && !quat_alg_elem_is_zero(&gen_rerand);
    }

    quat_alg_mul(&gen, &gen, &gen_rerand, (params->algebra));
    // in both cases, whether norm is prime or not prime,
    // gen is not divisible by any integer factor of the target norm
    // therefore the call below will yield an ideal of the correct norm
    quat_lideal_create(lideal, &gen, norm, &((params->order)->order), (params->algebra));
    assert(ibz_cmp(norm, &(lideal->norm)) == 0);

    ibz_finalize(&n_temp);
    quat_alg_elem_finalize(&gen);
    quat_alg_elem_finalize(&gen_rerand);
    ibz_finalize(&norm_d);
    ibz_finalize(&disc);
    return (found);
}

int quat_constant_random_O0_ideal_norm_1mod4(quat_left_ideal_t *lideal,
                                             const quat_represent_integer_params_t *params,
                                             const ibz_t *norm,
                                             const ibz_t *j,
                                             const ibz_t *x0,
                                             const ibz_t *y0)
{
    ibz_t x, y;
    ibz_t A, A_inv, u, v, sum, diff, two_j, inv_2j;
    ibz_t k1, k2, k3, u_plus_v, x0_plus_y0, N_minus_1;
    ibz_t two_N, two_N_minus_y;

    ibz_init(&x); ibz_init(&y);
    ibz_init(&A); ibz_init(&A_inv); ibz_init(&u); ibz_init(&v); 
    ibz_init(&sum); ibz_init(&diff); ibz_init(&two_j); ibz_init(&inv_2j);
    ibz_init(&k1); ibz_init(&k2); ibz_init(&k3); 
    ibz_init(&u_plus_v); ibz_init(&x0_plus_y0); ibz_init(&N_minus_1);
    ibz_init(&two_N); ibz_init(&two_N_minus_y);

    ibz_sub(&N_minus_1, norm, &ibz_const_one);
    ibz_rand_interval(&A, &ibz_const_one, &N_minus_1); //  A in [1, N-1]

    // A_inv = A^-1 mod N
    ibz_invmod(&A_inv, &A, norm);

    // Compute u = (A + A_inv)/2 mod N
    ibz_add(&sum, &A, &A_inv);
    if (ibz_is_odd(&sum)) {
        ibz_add(&sum, &sum, norm); // If (A + A_inv) is odd, then (A + A_inv + N) is even
    }
    ibz_div_2exp(&u, &sum, 1);     // Shift right by 1 bit (divide by 2)
    // (A + A_inv + N)/2 <= 1.5N
    // Use subtraction instead of modulo.
    if (ibz_cmp(&u, norm) >= 0) {
        ibz_sub(&u, &u, norm);
    }

    // Compute v = (A - A_inv)/(2*j) mod N
    ibz_add(&two_j, j, j);
    ibz_invmod(&inv_2j, &two_j, norm);
    ibz_sub(&diff, &A, &A_inv);
    ibz_mul(&v, &diff, &inv_2j);
    ibz_mod(&v, &v, norm); 

    // The original formula requires 4 large number multiplications:
    // x = x0*u - y0*v mod N
    // y = x0*v + y0*u mod N
    // Optimization: Gauss's complex multiplication technique, requires only 3 large number multiplications
    ibz_mul(&k1, x0, &u);                    // k1 = x0 * u
    ibz_mul(&k2, y0, &v);                    // k2 = y0 * v

    ibz_add(&u_plus_v, &u, &v);              // u + v
    ibz_add(&x0_plus_y0, x0, y0);            // x0 + y0
    ibz_mul(&k3, &x0_plus_y0, &u_plus_v);    // k3 = (x0 + y0)(u + v)

    // x = k1 - k2 mod N
    ibz_sub(&x, &k1, &k2);
    ibz_mod(&x, &x, norm);                  

    // y = k3 - k1 - k2 mod N
    ibz_sub(&y, &k3, &k1);
    ibz_sub(&y, &y, &k2);
    ibz_mod(&y, &y, norm);               

    // x, y in [0, 2*N], x is even and y is odd, make sure basis of I in O0 
    if (ibz_is_odd(&x)) { 
        ibz_add(&x, norm, &x); 
    }
    if (ibz_is_even(&y)) { 
        ibz_add(&y, norm, &y); 
    }

    // Generates left O0-ideals with norm N
    // Set lattices (reduces denominators)
    ibz_set(&lideal->lattice.denom, 2);
    ibz_add(&two_N, norm, norm);         // two_N = 2 * N
    ibz_sub(&two_N_minus_y, &two_N, &y); // two_N_minus_y = 2N - y
    ibz_mat_4x4_zero(&lideal->lattice.basis);
    // The basis matrix of lattices is
    //    alpha_1  alpha_2 alpha_3 alpha_4
    // 1 [ 1/2,       0,     0,      0    ]
    // i [  0,      1/2,     0,      0    ]
    // j [ x/2,     y/2,     N,      0    ]
    // k [(2N-y)/2, x/2,     0,      N    ]
    // column 0: 1/2 + (x/2)*j + ((2N-y)/2)*k -> [1, 0, x, 2N-y]^T
    ibz_set(&lideal->lattice.basis[0][0], 1);
    ibz_copy(&lideal->lattice.basis[2][0], &x);
    ibz_copy(&lideal->lattice.basis[3][0], &two_N_minus_y);
    // column 1: (1/2)*i + (y/2)*j + (x/2)*k -> [0, 1, y, x]^T
    ibz_set(&lideal->lattice.basis[1][1], 1);
    ibz_copy(&lideal->lattice.basis[2][1], &y);
    ibz_copy(&lideal->lattice.basis[3][1], &x);
    // column 2: N*j -> [0, 0, 2N, 0]^T
    ibz_copy(&lideal->lattice.basis[2][2], &two_N);
    // column 3: N*k -> [0, 0, 0, 2N]^T
    ibz_copy(&lideal->lattice.basis[3][3], &two_N);

    // Set order
    lideal->parent_order = &((params->order)->order);
    // Set norm
    ibz_copy(&lideal->norm, norm);

    ibz_finalize(&x); ibz_finalize(&y);
    ibz_finalize(&A); ibz_finalize(&A_inv); ibz_finalize(&u); ibz_finalize(&v); 
    ibz_finalize(&sum); ibz_finalize(&diff); ibz_finalize(&two_j); ibz_finalize(&inv_2j);
    ibz_finalize(&k1); ibz_finalize(&k2); ibz_finalize(&k3); 
    ibz_finalize(&u_plus_v); ibz_finalize(&x0_plus_y0); ibz_finalize(&N_minus_1);
    ibz_finalize(&two_N); ibz_finalize(&two_N_minus_y);

    return 1;
}

int 
quat_constant_random_O0_ideal_norm_3mod4(quat_left_ideal_t *lideal,
                                               const quat_represent_integer_params_t *params,
                                               const ibz_t *norm,
                                               const ibz_t *u0,
                                               const ibz_t *v0,
                                               const ibz_t *x0,
                                               const ibz_t *y0)
{
    ibz_t x, y;
    ibz_t k, u, v;
    ibz_t k1, k2, k3;
    ibz_t two_N, two_N_minus_y;

    ibz_init(&x); ibz_init(&y);
    ibz_init(&k); ibz_init(&u); ibz_init(&v);
    ibz_init(&k1); ibz_init(&k2); ibz_init(&k3);
    ibz_init(&two_N); ibz_init(&two_N_minus_y);

    // u + i*v = (u0 + i*v0)^k
    ibz_rand_interval(&k, &ibz_const_zero, norm); //  k in [0, N]
    ibz_complex_pow_mod(&u, &v, u0, v0, &k, norm); 

    // The original formula requires 4 large number multiplications:
    // x = x0*u - y0*v mod N
    // y = x0*v + y0*u mod N
    // Optimization: Gauss's complex multiplication technique, requires only 3 large number multiplications
    ibz_mul(&k1, x0, &u);   // k1 = x0 * u
    ibz_mul(&k2, y0, &v);   // k2 = y0 * v
    ibz_add(&x, &u, &v);
    ibz_add(&y, x0, y0);
    ibz_mul(&k3, &x, &y);   // k3 = (x0 + y0) * (u + v)

    // x = k1 - k2 mod N
    ibz_sub(&x, &k1, &k2);
    ibz_mod(&x, &x, norm);

    // y = k3 - k1 - k2 mod N
    ibz_sub(&y, &k3, &k1);
    ibz_sub(&y, &y, &k2);
    ibz_mod(&y, &y, norm);

    // x, y in [0, 2*N], x is even and y is odd, make sure basis of I in O0 
    if (ibz_is_odd(&x)) { 
        ibz_add(&x, norm, &x); 
    }
    if (ibz_is_even(&y)) { 
        ibz_add(&y, norm, &y); 
    }

    // Generates left O0-ideals with norm N
    // Set lattices (reduces denominators)
    ibz_set(&lideal->lattice.denom, 2);
    ibz_add(&two_N, norm, norm);         // two_N = 2 * N
    ibz_sub(&two_N_minus_y, &two_N, &y); // two_N_minus_y = 2N - y
    ibz_mat_4x4_zero(&lideal->lattice.basis);
    // The basis matrix of lattices is
    //    alpha_1  alpha_2 alpha_3 alpha_4
    // 1 [ 1/2,       0,     0,      0    ]
    // i [  0,      1/2,     0,      0    ]
    // j [ x/2,     y/2,     N,      0    ]
    // k [(2N-y)/2, x/2,     0,      N    ]
    // column 0: 1/2 + (x/2)*j + ((2N-y)/2)*k -> [1, 0, x, 2N-y]^T
    ibz_set(&lideal->lattice.basis[0][0], 1);
    ibz_copy(&lideal->lattice.basis[2][0], &x);
    ibz_copy(&lideal->lattice.basis[3][0], &two_N_minus_y);
    // column 1: (1/2)*i + (y/2)*j + (x/2)*k -> [0, 1, y, x]^T
    ibz_set(&lideal->lattice.basis[1][1], 1);
    ibz_copy(&lideal->lattice.basis[2][1], &y);
    ibz_copy(&lideal->lattice.basis[3][1], &x);
    // column 2: N*j -> [0, 0, 2N, 0]^T
    ibz_copy(&lideal->lattice.basis[2][2], &two_N);
    // column 3: N*k -> [0, 0, 0, 2N]^T
    ibz_copy(&lideal->lattice.basis[3][3], &two_N);

    // Set order
    lideal->parent_order = &((params->order)->order);
    // Set norm
    ibz_copy(&lideal->norm, norm);

    ibz_finalize(&x); ibz_finalize(&y);
    ibz_finalize(&k); ibz_finalize(&u); ibz_finalize(&v);
    ibz_finalize(&k1); ibz_finalize(&k2); ibz_finalize(&k3);
    ibz_finalize(&two_N); ibz_finalize(&two_N_minus_y);

    return 1;
}

int quat_constant_random_O0_ideal_norm_universal(quat_left_ideal_t *lideal,
                                                   const quat_represent_integer_params_t *params,
                                                   const ibz_t *norm,
                                                   const ibz_t *x0,
                                                   const ibz_t *y0,
                                                   const ibz_t *minus_y0,
                                                   const ibz_t *x0_plus_y0)
{
    ibz_t x, y;
    ibz_t t, t2, den, inv_den;
    ibz_t A, B, A_minus_B;
    ibz_t k1, k2, k3, num_x, num_y;
    ibz_t two_N, two_N_minus_y;
    int is_fallback = 0;

    ibz_init(&x); ibz_init(&y);
    ibz_init(&t); ibz_init(&t2); ibz_init(&den); ibz_init(&inv_den);
    ibz_init(&A); ibz_init(&B); ibz_init(&A_minus_B);
    ibz_init(&k1); ibz_init(&k2); ibz_init(&k3); ibz_init(&num_x); ibz_init(&num_y);
    ibz_init(&two_N); ibz_init(&two_N_minus_y);

    ibz_rand_interval(&t, &ibz_const_zero, norm); //  t in [0, N]

    // Compute t2 = t^2 mod N, den = (t^2 + 1) mod N, inv_den = 1/(t^2 + 1) mod N
    if (ibz_cmp(&t, norm) == 0) {
        is_fallback = 1;
    } 
    else {
        // t2 = t^2 mod N, t2 in [0, N-1]
        ibz_mul(&t2, &t, &t);
        ibz_mod(&t2, &t2, norm);

        // den = t2 + 1 mod N
        ibz_add(&den, &t2, &ibz_const_one);

        // check if den == N (Can only happen when N = 1 mod 4)
        if (ibz_cmp(&den, norm) == 0) {
            is_fallback = 1;
        } else {
            // inv_den = 1/den mod N
            if (!ibz_invmod(&inv_den, &den, norm)) {
                is_fallback = 1; 
            }
        }
    }

    // Compute x = (x0*(t^2-1) + 2*y0*t)/(t^2+1) mod N
    //         y = (y0*(1-t^2) + 2*x0*t)/(t^2+1) mod N
    if (is_fallback) {
        // Can only happen when t = N or t^2 + 1 = 0 mod N
        // (x, y) is the point symmetric to (x0, y0) about the x-axis.
        ibz_copy(&x, x0);
        ibz_copy(&y, minus_y0);
    } 
    else {
        // A = t^2 - 1
        ibz_sub(&A, &t2, &ibz_const_one);
        // B = 2*t
        ibz_add(&B, &t, &t);
        // A_minus_B = A - B
        ibz_sub(&A_minus_B, &A, &B);

        // The original formula requires 4 large number multiplications:
        // num_x = x0 * A + y0 * B 
        // num_y = x0 * B - y0 * A 
        // Optimization: Gauss's complex multiplication technique, requires only 3 large number multiplications
        ibz_mul(&k1, x0, &A);                   // k1 = x0 * A
        ibz_mul(&k2, y0, &B);                   // k2 = y0 * B
        ibz_mul(&k3, x0_plus_y0, &A_minus_B);   // k3 = (x0 + y0) * (A - B)

        // num_x = k1 + k2  (= x0*A + y0*B)
        ibz_add(&num_x, &k1, &k2);

        // num_y = k1 - k2 - k3 (= x0*B - y0*A)
        ibz_sub(&num_y, &k1, &k2);
        ibz_sub(&num_y, &num_y, &k3);

        // x = (num_x * inv_den) mod N
        ibz_mul(&x, &num_x, &inv_den);
        ibz_mod(&x, &x, norm);

        // y = (num_y * inv_den) mod N
        ibz_mul(&y, &num_y, &inv_den);
        ibz_mod(&y, &y, norm);
    }
    
    // x, y in [0, 2*N], x is even and y is odd, make sure basis of I in O0 
    if (ibz_is_odd(&x)) { 
        ibz_add(&x, norm, &x); 
    }
    if (ibz_is_even(&y)) { 
        ibz_add(&y, norm, &y); 
    }

    // Generates left O0-ideals with norm N
    // Set lattices (reduces denominators)
    ibz_set(&lideal->lattice.denom, 2);
    ibz_add(&two_N, norm, norm);         // two_N = 2 * N
    ibz_sub(&two_N_minus_y, &two_N, &y); // two_N_minus_y = 2N - y
    ibz_mat_4x4_zero(&lideal->lattice.basis);
    // The basis matrix of lattices is
    //    alpha_1  alpha_2 alpha_3 alpha_4
    // 1 [ 1/2,       0,     0,      0    ]
    // i [  0,      1/2,     0,      0    ]
    // j [ x/2,     y/2,     N,      0    ]
    // k [(2N-y)/2, x/2,     0,      N    ]
    // column 0: 1/2 + (x/2)*j + ((2N-y)/2)*k -> [1, 0, x, 2N-y]^T
    ibz_set(&lideal->lattice.basis[0][0], 1);
    ibz_copy(&lideal->lattice.basis[2][0], &x);
    ibz_copy(&lideal->lattice.basis[3][0], &two_N_minus_y);
    // column 1: (1/2)*i + (y/2)*j + (x/2)*k -> [0, 1, y, x]^T
    ibz_set(&lideal->lattice.basis[1][1], 1);
    ibz_copy(&lideal->lattice.basis[2][1], &y);
    ibz_copy(&lideal->lattice.basis[3][1], &x);
    // column 2: N*j -> [0, 0, 2N, 0]^T
    ibz_copy(&lideal->lattice.basis[2][2], &two_N);
    // column 3: N*k -> [0, 0, 0, 2N]^T
    ibz_copy(&lideal->lattice.basis[3][3], &two_N);

    // Set order
    lideal->parent_order = &((params->order)->order);
    // Set norm
    ibz_copy(&lideal->norm, norm);

    ibz_finalize(&x); ibz_finalize(&y);
    ibz_finalize(&t); ibz_finalize(&t2); ibz_finalize(&den); ibz_finalize(&inv_den);
    ibz_finalize(&A); ibz_finalize(&B); ibz_finalize(&A_minus_B);
    ibz_finalize(&k1); ibz_finalize(&k2); ibz_finalize(&k3); ibz_finalize(&num_x); ibz_finalize(&num_y);
    ibz_finalize(&two_N); ibz_finalize(&two_N_minus_y);

    return 1;
}

int quat_generate_random_O0_ideal_norm_non_prime(quat_left_ideal_t *lideal,
                                                 const quat_represent_integer_params_t *params,
                                                 const ibz_t *norm,
                                                 const ibz_t *primes,     
                                                 const ibz_t *powers,    
                                                 int num_factors)
{
    if (num_factors <= 0) return 0;

    ibz_t c, xi, yi, m, target, inv_y, inv_2, tmp, two, q_minus_1, M_i, inv_Mi, rem, two_N, two_N_minus_y;
    ibz_t x, y;
    ibz_init(&c); ibz_init(&xi); ibz_init(&yi); ibz_init(&m); ibz_init(&target); 
    ibz_init(&inv_y); ibz_init(&inv_2); ibz_init(&tmp); ibz_init(&two); 
    ibz_init(&q_minus_1); ibz_init(&M_i); ibz_init(&inv_Mi); ibz_init(&rem);
    ibz_init(&x); ibz_init(&y);
    ibz_init(&two_N); ibz_init(&two_N_minus_y);

    ibz_copy(&two, &ibz_const_two);
    ibz_copy(&x, &ibz_const_zero);
    ibz_copy(&y, &ibz_const_zero);

    // Compute c = -p^{-1} mod N
    ibz_invmod(&c, &(params->algebra->p), norm);
    ibz_sub(&c, norm, &c);

    // Compute x^2 + y^2 = c mod N
    for (int i = 0; i < num_factors; i++) {
        ibz_sub(&q_minus_1, &primes[i], &ibz_const_one); 
        
        // Tonelli-Shanks: compute x_i^2 + y_i^2 = c mod q_i
        while (1) {
            ibz_rand_interval(&xi, &ibz_const_zero, &q_minus_1); // x_i in [0, q_i - 1]
            ibz_mul(&tmp, &xi, &xi);
            ibz_sub(&tmp, &c, &tmp);          
            ibz_mod(&tmp, &tmp, &primes[i]); 
            
            if (ibz_is_zero(&tmp) && ibz_cmp(&primes[i], &powers[i]) != 0) continue; 
            if (ibz_legendre(&tmp, &primes[i]) == 1) {
                ibz_sqrt_mod_p(&yi, &tmp, &primes[i]); 
                break;
            }
        }

        // Hensel boost: compute x_i^2 + y_i^2 = c mod q_i^{e_i} 
        // if e_i = 1, skip
        if (ibz_cmp(&primes[i], &powers[i]) != 0) {
            ibz_copy(&m, &primes[i]);
            ibz_mul(&target, &xi, &xi);
            ibz_sub(&target, &c, &target); 
            ibz_mod(&target, &target, &powers[i]); // Target = c - x_i^2 mod q_i^{e_i}

            while (ibz_cmp(&m, &powers[i]) < 0) {
                ibz_mul(&m, &m, &m); 
                if (ibz_cmp(&m, &powers[i]) > 0) ibz_copy(&m, &powers[i]); 

                // y_new = (y_old + target * y_old^(-1)) * 2^(-1) mod q_i
                // After each iteration, the exponent doubles
                ibz_invmod(&inv_y, &yi, &m);
                ibz_mul(&tmp, &target, &inv_y);
                ibz_add(&yi, &yi, &tmp);

                ibz_invmod(&inv_2, &two, &m);
                ibz_mul(&yi, &yi, &inv_2);
                ibz_mod(&yi, &yi, &m);
            }
        }

        // CRT: compute x = x_1 * M_1 * M_1^{-1} + ... + x_n * M_n * M_n^{-1}
        //              y = y_1 * M_1 * M_1^{-1} + ... + y_n * M_n * M_n^{-1}
        // x^2 + y^2 = c mod N
        ibz_div(&M_i, &rem, norm, &powers[i]); // M_i = N / q_i^{e_i}
        ibz_invmod(&inv_Mi, &M_i, &powers[i]); // M_i^-1 mod q_i^{e_i}

        // out_x += x_i * M_i * inv_Mi
        ibz_mul(&tmp, &xi, &M_i);
        ibz_mul(&tmp, &tmp, &inv_Mi);
        ibz_add(&x, &x, &tmp);

        // out_y += y_i * M_i * inv_Mi
        ibz_mul(&tmp, &yi, &M_i);
        ibz_mul(&tmp, &tmp, &inv_Mi);
        ibz_add(&y, &y, &tmp);
    }
    // Global modulo N
    ibz_mod(&x, &x, norm);
    ibz_mod(&y, &y, norm);

    // x, y in [0, 2*N], x is even and y is odd, make sure basis of I in O0 
    if (ibz_is_odd(&x)) { 
        ibz_add(&x, norm, &x); 
    }
    if (ibz_is_even(&y)) { 
        ibz_add(&y, norm, &y); 
    }

    // Generates left O0-ideals with norm N
    // Set lattices (reduces denominators)
    ibz_set(&lideal->lattice.denom, 2);
    ibz_add(&two_N, norm, norm);         // two_N = 2 * N
    ibz_sub(&two_N_minus_y, &two_N, &y); // two_N_minus_y = 2N - y
    ibz_mat_4x4_zero(&lideal->lattice.basis);
    // The basis matrix of lattices is
    //    alpha_1  alpha_2 alpha_3 alpha_4
    // 1 [ 1/2,       0,     0,      0    ]
    // i [  0,      1/2,     0,      0    ]
    // j [ x/2,     y/2,     N,      0    ]
    // k [(2N-y)/2, x/2,     0,      N    ]
    // column 0: 1/2 + (x/2)*j + ((2N-y)/2)*k -> [1, 0, x, 2N-y]^T
    ibz_set(&lideal->lattice.basis[0][0], 1);
    ibz_copy(&lideal->lattice.basis[2][0], &x);
    ibz_copy(&lideal->lattice.basis[3][0], &two_N_minus_y);
    // column 1: (1/2)*i + (y/2)*j + (x/2)*k -> [0, 1, y, x]^T
    ibz_set(&lideal->lattice.basis[1][1], 1);
    ibz_copy(&lideal->lattice.basis[2][1], &y);
    ibz_copy(&lideal->lattice.basis[3][1], &x);
    // column 2: N*j -> [0, 0, 2N, 0]^T
    ibz_copy(&lideal->lattice.basis[2][2], &two_N);
    // column 3: N*k -> [0, 0, 0, 2N]^T
    ibz_copy(&lideal->lattice.basis[3][3], &two_N);

    // Set order
    lideal->parent_order = &((params->order)->order);
    // Set norm
    ibz_copy(&lideal->norm, norm);

    ibz_finalize(&c); ibz_finalize(&xi); ibz_finalize(&yi); ibz_finalize(&m); 
    ibz_finalize(&target); ibz_finalize(&inv_y); ibz_finalize(&inv_2); 
    ibz_finalize(&tmp); ibz_finalize(&two); ibz_finalize(&q_minus_1);
    ibz_finalize(&M_i); ibz_finalize(&inv_Mi); ibz_finalize(&rem);
    ibz_finalize(&x); ibz_finalize(&y);
    ibz_finalize(&two_N); ibz_finalize(&two_N_minus_y);

    return 1;
}

void
quat_change_to_O0_basis(ibz_vec_4_t *vec, const quat_alg_elem_t *el)
{
    ibz_t tmp;
    ibz_init(&tmp);
    ibz_copy(&(*vec)[2], &el->coord[2]);
    ibz_add(&(*vec)[2], &(*vec)[2], &(*vec)[2]); // double (not optimal if el->denom is even...)
    ibz_copy(&(*vec)[3], &el->coord[3]);         // double (not optimal if el->denom is even...)
    ibz_add(&(*vec)[3], &(*vec)[3], &(*vec)[3]);
    ibz_sub(&(*vec)[0], &el->coord[0], &el->coord[3]);
    ibz_sub(&(*vec)[1], &el->coord[1], &el->coord[2]);

    assert(ibz_divides(&(*vec)[0], &el->denom));
    assert(ibz_divides(&(*vec)[1], &el->denom));
    assert(ibz_divides(&(*vec)[2], &el->denom));
    assert(ibz_divides(&(*vec)[3], &el->denom));

    ibz_div(&(*vec)[0], &tmp, &(*vec)[0], &el->denom);
    ibz_div(&(*vec)[1], &tmp, &(*vec)[1], &el->denom);
    ibz_div(&(*vec)[2], &tmp, &(*vec)[2], &el->denom);
    ibz_div(&(*vec)[3], &tmp, &(*vec)[3], &el->denom);

    ibz_finalize(&tmp);
}

