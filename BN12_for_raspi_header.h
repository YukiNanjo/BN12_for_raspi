#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#define x_bit 114
#define x_bit_for_opt_ate 116
int X_bit_binary[x_bit+1];
int X_bit_binary_for_opt_ate[x_bit_for_opt_ate+1];

mpz_t X;
mpz_t prime,EFp_order,trace_t;
mpz_t EFp_total,EFp12_total;
mpz_t curve_b;
mpz_t final_exp;

struct Fp{
    mpz_t x0;
};

struct Fp2{
    struct Fp x0;
    struct Fp x1;
};

struct Fp4{
    struct Fp2 x0;
    struct Fp2 x1;
};

struct Fp12{
    struct Fp4 x0;
    struct Fp4 x1;
    struct Fp4 x2;
};

struct EFp{
    struct Fp x;
    struct Fp y;
    int infinity;
};

struct EFp2{
    struct Fp2 x;
    struct Fp2 y;
    int infinity;
};

struct EFp4{
    struct Fp4 x;
    struct Fp4 y;
    int infinity;
};

struct EFp12{
    struct Fp12 x;
    struct Fp12 y;
    int infinity;
};

struct Fp2 Fp2_basis,Fp2_basis_inv;
struct Fp epsilon1,epsilon2;
enum state{
    f_p1,f_p2,f_p3,f_p4,f_p5,f_p6,f_p7,f_p8,f_p9,f_p10,f_p11,f_p12
};
struct Fp2 frobenius_constant[12][6];
struct Fp2 skew_frobenius_constant[12][2];

struct timeval t0,t1;
float OPT_ATE_MILLER,X_ATE_MILLER;
float PLAIN_FINAL_EXP,OPT_FINAL_EXP1,OPT_FINAL_EXP2;
float OPT_ATE_TOTAL,X_ATE_TOTAL;
float PLAIN_G1_SCM,G1_SCM_2SPLIT;
float PLAIN_G2_SCM,G2_SCM_2SPLIT,G2_SCM_4SPLIT;
float PLAIN_G3_EXP,G3_EXP_2SPLIT,G3_EXP_4SPLIT;

/*============================================================================*/
/* field                                                                      */
/*============================================================================*/
//Fp
void Fp_init(struct Fp *P);
void Fp_set(struct Fp *P,struct Fp *A);
void Fp_set_ui(struct Fp *P,unsigned long int a);
void Fp_set_mpz(struct Fp *P,mpz_t a);
void Fp_set_neg(struct Fp *P,struct Fp *A);
void Fp_random(struct Fp *P,gmp_randstate_t state);
void Fp_clear(struct Fp *P);
void Fp_printf(struct Fp *P,char *name);
void Fp_mul(struct Fp *ANS,struct Fp *A,struct Fp *B);
void Fp_mul_ui(struct Fp *ANS,struct Fp *A,unsigned long int a);
void Fp_mul_mpz(struct Fp *ANS,struct Fp *A,mpz_t a);
void Fp_add(struct Fp *ANS,struct Fp *A,struct Fp *B);
void Fp_add_ui(struct Fp *ANS,struct Fp *A,unsigned long int a);
void Fp_add_mpz(struct Fp *ANS,struct Fp *A,mpz_t a);
void Fp_sub(struct Fp *ANS,struct Fp *A,struct Fp *B);
void Fp_sub_ui(struct Fp *ANS,struct Fp *A,unsigned long int a);
void Fp_sub_mpz(struct Fp *ANS,struct Fp *A,mpz_t a);
void Fp_inv(struct Fp *ANS,struct Fp *A);
int Fp_legendre(struct Fp *A);
int Fp_isCNR(struct Fp *A);
void Fp_sqrt(struct Fp *ANS,struct Fp *A);
void Fp_pow(struct Fp *ANS,struct Fp *A,mpz_t a);
int Fp_cmp(struct Fp *A,struct Fp *B);
int Fp_cmp_ui(struct Fp *A,unsigned long int a);
int Fp_cmp_mpz(struct Fp *A,mpz_t a);
int Fp_cmp_zero(struct Fp *A);
int Fp_cmp_one(struct Fp *A);
//Fp2
void Fp2_init(struct Fp2 *P);
void Fp2_set(struct Fp2 *P,struct Fp2 *A);
void Fp2_set_ui(struct Fp2 *P,unsigned long int a);
void Fp2_set_mpz(struct Fp2 *P,mpz_t a);
void Fp2_set_neg(struct Fp2 *P,struct Fp2 *A);
void Fp2_random(struct Fp2 *P,gmp_randstate_t state);
void Fp2_clear(struct Fp2 *P);
void Fp2_printf(struct Fp2 *P,char *name);
void Fp2_mul(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B);
void Fp2_mul_ui(struct Fp2 *ANS,struct Fp2 *A,unsigned long int a);
void Fp2_mul_mpz(struct Fp2 *ANS,struct Fp2 *A,mpz_t a);
void Fp2_mul_basis(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_sqr(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_inv_basis(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_add(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B);
void Fp2_add_ui(struct Fp2 *ANS,struct Fp2 *A,unsigned long int a);
void Fp2_add_mpz(struct Fp2 *ANS,struct Fp2 *A,mpz_t a);
void Fp2_sub(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B);
void Fp2_sub_ui(struct Fp2 *ANS,struct Fp2 *A,unsigned long int a);
void Fp2_sub_mpz(struct Fp2 *ANS,struct Fp2 *A,mpz_t a);
void Fp2_inv(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_inv_map(struct Fp2 *ANS,struct Fp2 *A);
int Fp2_legendre(struct Fp2 *A);
void Fp2_sqrt(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_pow(struct Fp2 *ANS,struct Fp2 *A,mpz_t a);
int Fp2_cmp(struct Fp2 *A,struct Fp2 *B);
int Fp2_cmp_ui(struct Fp2 *A,unsigned long int a);
int Fp2_cmp_mpz(struct Fp2 *A,mpz_t a);
int Fp2_cmp_zero(struct Fp2 *A);
int Fp2_cmp_one(struct Fp2 *A);
//Fp4
void Fp4_init(struct Fp4 *P);
void Fp4_set(struct Fp4 *P,struct Fp4 *A);
void Fp4_set_ui(struct Fp4 *P,unsigned long int a);
void Fp4_set_mpz(struct Fp4 *P,mpz_t a);
void Fp4_set_neg(struct Fp4 *P,struct Fp4 *A);
void Fp4_random(struct Fp4 *P,gmp_randstate_t state);
void Fp4_clear(struct Fp4 *P);
void Fp4_printf(struct Fp4 *P,char *name);
void Fp4_mul(struct Fp4 *ANS,struct Fp4 *A,struct Fp4 *B);
void Fp4_mul_ui(struct Fp4 *ANS,struct Fp4 *A,unsigned long int a);
void Fp4_mul_mpz(struct Fp4 *ANS,struct Fp4 *A,mpz_t a);
void Fp4_mul_basis(struct Fp4 *ANS,struct Fp4 *A);
void Fp4_sqr(struct Fp4 *ANS,struct Fp4 *A);
void Fp4_add(struct Fp4 *ANS,struct Fp4 *A,struct Fp4 *B);
void Fp4_add_ui(struct Fp4 *ANS,struct Fp4 *A,unsigned long int a);
void Fp4_add_mpz(struct Fp4 *ANS,struct Fp4 *A,mpz_t a);
void Fp4_sub(struct Fp4 *ANS,struct Fp4 *A,struct Fp4 *B);
void Fp4_sub_ui(struct Fp4 *ANS,struct Fp4 *A,unsigned long int a);
void Fp4_sub_mpz(struct Fp4 *ANS,struct Fp4 *A,mpz_t a);
void Fp4_inv(struct Fp4 *ANS,struct Fp4 *A);
void Fp4_inv_map(struct Fp4 *ANS,struct Fp4 *A);
int Fp4_legendre(struct Fp4 *A);
void Fp4_sqrt(struct Fp4 *ANS,struct Fp4 *A);
void Fp4_pow(struct Fp4 *ANS,struct Fp4 *A,mpz_t a);
int Fp4_cmp(struct Fp4 *A,struct Fp4 *B);
int Fp4_cmp_ui(struct Fp4 *A,unsigned long int a);
int Fp4_cmp_mpz(struct Fp4 *A,mpz_t a);
int Fp4_cmp_zero(struct Fp4 *A);
int Fp4_cmp_one(struct Fp4 *A);
//Fp12
void Fp12_init(struct Fp12 *P);
void Fp12_set(struct Fp12 *P,struct Fp12 *A);
void Fp12_set_ui(struct Fp12 *P,unsigned long int a);
void Fp12_set_mpz(struct Fp12 *P,mpz_t a);
void Fp12_set_neg(struct Fp12 *P,struct Fp12 *A);
void Fp12_random(struct Fp12 *P,gmp_randstate_t state);
void Fp12_clear(struct Fp12 *P);
void Fp12_printf(struct Fp12 *P,char *name);
void Fp12_mul(struct Fp12 *ANS,struct Fp12 *A,struct Fp12 *B);
void Fp12_mul_ui(struct Fp12 *ANS,struct Fp12 *A,unsigned long int a);
void Fp12_mul_mpz(struct Fp12 *ANS,struct Fp12 *A,mpz_t a);
void Fp12_sqr(struct Fp12 *ANS,struct Fp12 *A);
void Fp12_add(struct Fp12 *ANS,struct Fp12 *A,struct Fp12 *B);
void Fp12_add_ui(struct Fp12 *ANS,struct Fp12 *A,unsigned long int a);
void Fp12_add_mpz(struct Fp12 *ANS,struct Fp12 *A,mpz_t a);
void Fp12_sub(struct Fp12 *ANS,struct Fp12 *A,struct Fp12 *B);
void Fp12_sub_ui(struct Fp12 *ANS,struct Fp12 *A,unsigned long int a);
void Fp12_sub_mpz(struct Fp12 *ANS,struct Fp12 *A,mpz_t a);
void Fp12_inv(struct Fp12 *ANS,struct Fp12 *A);
void Fp12_inv_map1(struct Fp12 *ANS,struct Fp12 *A);
void Fp12_inv_map2(struct Fp12 *ANS,struct Fp12 *A);
int Fp12_legendre(struct Fp12 *A);
void Fp12_sqrt(struct Fp12 *ANS,struct Fp12 *A);
void Fp12_pow(struct Fp12 *ANS,struct Fp12 *A,mpz_t a);
int Fp12_cmp(struct Fp12 *A,struct Fp12 *B);
int Fp12_cmp_ui(struct Fp12 *A,unsigned long int a);
int Fp12_cmp_mpz(struct Fp12 *A,mpz_t a);
int Fp12_cmp_zero(struct Fp12 *A);
int Fp12_cmp_one(struct Fp12 *A);
void Fp12_frobenius_map_p1(struct Fp12 *ANS,struct Fp12 *A);
void Fp12_frobenius_map_p2(struct Fp12 *ANS,struct Fp12 *A);
void Fp12_frobenius_map_p3(struct Fp12 *ANS,struct Fp12 *A);
void Fp12_frobenius_map_p4(struct Fp12 *ANS,struct Fp12 *A);
void Fp12_frobenius_map_p6(struct Fp12 *ANS,struct Fp12 *A);
void Fp12_frobenius_map_p8(struct Fp12 *ANS,struct Fp12 *A);
void Fp12_frobenius_map_p10(struct Fp12 *ANS,struct Fp12 *A);
/*============================================================================*/
/* Elliptic curve                                                             */
/*============================================================================*/
//EFp
void EFp_init(struct EFp *P);
void EFp_set(struct EFp *P,struct EFp *A);
void EFp_set_ui(struct EFp *P,unsigned long int a);
void EFp_set_mpz(struct EFp *P,mpz_t a);
void EFp_set_neg(struct EFp *P,struct EFp *A);
void EFp_clear(struct EFp *P);
void EFp_printf(struct EFp *P,char *name);
void EFp_rational_point(struct EFp *P);
void EFp_ECD(struct EFp *ANS,struct EFp *P);
void EFp_ECA(struct EFp *ANS,struct EFp *P1,struct EFp *P2);
void EFp_SCM(struct EFp *ANS,struct EFp *P,mpz_t R);
void EFp_skew_frobenius_map_p2(struct EFp *ANS,struct EFp *A);
//EFp2
void EFp2_init(struct EFp2 *P);
void EFp2_set(struct EFp2 *P,struct EFp2 *A);
void EFp2_set_ui(struct EFp2 *P,unsigned long int a);
void EFp2_set_mpz(struct EFp2 *P,mpz_t a);
void EFp2_set_neg(struct EFp2 *ANS,struct EFp2 *P);
void EFp2_clear(struct EFp2 *P);
void EFp2_printf(struct EFp2 *P,char *name);
void EFp2_rational_point(struct EFp2 *P);
void EFp2_ECD(struct EFp2 *ANS,struct EFp2 *P);
void EFp2_ECA(struct EFp2 *ANS,struct EFp2 *P1,struct EFp2 *P2);
void EFp2_SCM(struct EFp2 *ANS,struct EFp2 *P,mpz_t R);
void EFp2_skew_frobenius_map_p1(struct EFp2 *ANS,struct EFp2 *A);
void EFp2_skew_frobenius_map_p2(struct EFp2 *ANS,struct EFp2 *A);
void EFp2_skew_frobenius_map_p3(struct EFp2 *ANS,struct EFp2 *A);
void EFp2_skew_frobenius_map_p10(struct EFp2 *ANS,struct EFp2 *A);
//EFp4
void EFp4_init(struct EFp4 *P);
void EFp4_set(struct EFp4 *P,struct EFp4 *A);
void EFp4_set_ui(struct EFp4 *P,unsigned long int a);
void EFp4_set_mpz(struct EFp4 *P,mpz_t a);
void EFp4_set_neg(struct EFp4 *P,struct EFp4 *A);
void EFp4_clear(struct EFp4 *P);
void EFp4_printf(struct EFp4 *P,char *name);
void EFp4_rational_point(struct EFp4 *P);
void EFp4_ECD(struct EFp4 *ANS,struct EFp4 *P);
void EFp4_ECA(struct EFp4 *ANS,struct EFp4 *P1,struct EFp4 *P2);
void EFp4_SCM(struct EFp4 *ANS,struct EFp4 *P,mpz_t R);
//EFp12
void EFp12_init(struct EFp12 *P);
void EFp12_set(struct EFp12 *P,struct EFp12 *A);
void EFp12_set_ui(struct EFp12 *P,unsigned long int a);
void EFp12_set_mpz(struct EFp12 *P,mpz_t a);
void EFp12_set_neg(struct EFp12 *P,struct EFp12 *A);
void EFp12_clear(struct EFp12 *P);
void EFp12_printf(struct EFp12 *P,char *name);
void EFp12_rational_point(struct EFp12 *P);
void EFp12_generate_G1(struct EFp12 *P);
void EFp12_generate_G2(struct EFp12 *Q);
void EFp12_ECD(struct EFp12 *ANS,struct EFp12 *P);
void EFp12_ECA(struct EFp12 *ANS,struct EFp12 *P1,struct EFp12 *P2);
void EFp12_SCM(struct EFp12 *ANS,struct EFp12 *P,mpz_t R);
/*============================================================================*/
/* pairing functions                                                          */
/*============================================================================*/
//G1 SCM
void EFp12_plain_G1_SCM(struct EFp12 *ANS,struct EFp12 *P,mpz_t S);
void EFp12_G1_SCM(struct EFp12 *ANS,struct EFp12 *P,mpz_t S);
//G2 SCM
void EFp12_plain_G2_SCM(struct EFp12 *ANS,struct EFp12 *Q,mpz_t S);
void EFp12_2split_G2_SCM(struct EFp12 *ANS,struct EFp12 *P,mpz_t S);
void EFp12_4split_G2_SCM(struct EFp12 *ANS,struct EFp12 *Q,mpz_t S);
//G3 EXP
void Fp12_plainG3_exp(struct Fp12 *ANS,struct Fp12 *A,mpz_t S);
void Fp12_2split_G3_exp(struct Fp12 *ANS,struct Fp12 *A,mpz_t S);
void Fp12_4split_G3_exp(struct Fp12 *ANS,struct Fp12 *A,mpz_t S);
//twist
void EFp12_to_EFp2(struct EFp2 *ANS,struct EFp12 *P);
void EFp2_to_EFp12(struct EFp12 *ANS,struct EFp2 *P);
//final exp
void EFp12_pairing_final_exp(struct Fp12 *ANS,struct Fp12 *A);              //normal
void EFp12_pairing_final_exp_optimal1(struct Fp12 *ANS,struct Fp12 *A);         //which I found
void EFp12_pairing_final_exp_optimal2(struct Fp12 *ANS,struct Fp12 *A);         //which Alamin san gave
//7-sparse
void Fp12_7_sparse_mul(struct Fp12 *ANS,struct Fp12 *A,struct Fp12 *B);
void ff_ltt_for_ate(struct Fp12 *f,struct EFp2 *T,struct EFp *P);
void f_ltq_for_ate(struct Fp12 *f,struct EFp2 *T,struct EFp2 *Q,struct EFp *P);
//miller
void miller_for_opt_ate(struct Fp12 *ANS,struct EFp12 *Q,struct EFp12 *P);
void miller_for_x_ate(struct Fp12 *ANS,struct EFp12 *Q,struct EFp12 *P);
//pairing
void Opt_ate_pairing(struct Fp12 *ANS,struct EFp12 *Q,struct EFp12 *P);
void X_ate_pairing(struct Fp12 *ANS,struct EFp12 *Q,struct EFp12 *P);
//init
void init_parameters();
//clear
void clear_parameters();
//set
void set_parameters();
void generate_X();
int generate_prime();
int generate_order();
void generate_trace();
void generate_basis();
void get_epsilon();
void get_scalar_of_final_exp();
void weil();
//print
void print_parameters();
//time
float timedifference_msec(struct timeval t0, struct timeval t1);
float timedifference_usec(struct timeval t0, struct timeval t1);
//test
void test();
void test_G1_SCM();
void test_G2_SCM();
void test_G3_pow();
void test_tate_pairing();
void test_plain_ate_pairing();
void test_opt_ate_pairing();
void test_x_ate_pairing();
void test_for_RaspberryPi();

