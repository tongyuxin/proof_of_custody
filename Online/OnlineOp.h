/*
Copyright (c) 2017, The University of Bristol, Senate House, Tyndall Avenue, Bristol, BS8 1TH, United Kingdom.
Copyright (c) 2019, COSIC-KU Leuven, Kasteelpark Arenberg 10, bus 2452, B-3001 Leuven-Heverlee, Belgium.

All rights reserved
*/
#ifndef _OnlineOP
#define _OnlineOP

#include "Offline/offline_data.h"
#include "Online/Machine.h"
#include "System/Player.h"
#include "Processor/Processor.h"
#include "LSSS/PRSS.h"

extern vector<sacrificed_data> SacrificeD;

enum
{
    TRIPLE = 0x50,
    BIT = 0x51,
    SQUARE = 0x52,
    INPUT_MASK = 0x53
};

class UsedTuples
{
public:
    unsigned int UsedTriples = 0;
    unsigned int UsedSquares = 0;
    unsigned int UsedBit = 0;
    unsigned int UsedInputMask = 0;
};

class Complex
{
public:
    Share real;
    Share imag;
    Complex() {}
    Complex(Share &_real, Share &_imag) : real(_real), imag(_imag) {}
    void setValue(Share &_real, Share &_imag)
    {
        real = _real;
        imag = _imag;
    }
};

class Complex_Plain
{
public:
    gfp real;
    gfp imag;
    Complex_Plain() {}
    Complex_Plain(gfp &_real, gfp &_imag) : real(_real), imag(_imag) {}
    void setValue(gfp &_real, gfp &_imag)
    {
        real = _real;
        imag = _imag;
    }
};

class OnlineOp
{
public:
    UsedTuples UT;
    Processor &Proc;
    int online_num;
    Player &P;
    offline_control_data &OCD;
    Machine &machine;
    PRSS prss;
    explicit OnlineOp(Processor &Proc_, int online_num_, Player &P_,
                      offline_control_data &OCD_, Machine &machine_)
        : Proc(Proc_), online_num(online_num_), P(P_), OCD(OCD_), machine(machine_)
    {
        prss = PRSS(P);
    }
    int verbose = 0;

    void getTuples(vector<Share> &sp, int opcode);

    /*Share ops*/
    // c = a + b (b is share)
    void add(Share &c, const Share &a, const Share &b);
    // c = a + b (b is plain)
    void add_plain(Share &c, const Share &a, const gfp &b);
    // c = a - b (b is share)
    void sub(Share &c, const Share &a, const Share &b);
    // c = a * b (b is plain)
    void mul_plain(Share &c, const Share &a, const gfp &b);
    // c = a * b (b is share)
    void mul(Share &c, const Share &a, const Share &b);
    // aa = a^2
    void sqr(Share &aa, const Share &a);
    //ia = a^{-1} mod q
    void inv(Share &ia, const Share &a);
    // c = a * b^{-1} mod q
    void div(Share &c, const Share &a, const Share &b);

    /*Complex ops*/
    // c = a + b (b is shared complex)
    void add(Complex &c, const Complex &a, const Complex &b);
    // c = a + b (b is plain complex)
    void add_plain(Complex &c, const Complex &a, const Complex_Plain &b);
    // c = a - b (b is shared complex)
    void sub(Complex &c, const Complex &a, const Complex &b);
    // c = a * b (b is plain complex)
    void mul_plain(Complex &c, const Complex &a, const Complex_Plain &b);
    // c = a * b (b is shared complex)
    void mul(Complex &c, const Complex &a, const Complex &b);
    // aa = a^2
    void sqr(Complex &aa, const Complex &a);
    //ia = a^{-1} mod (q,x^2+1)
    void inv(Complex &ia, const Complex a);
    // c = a * b^{-1} mod (q, x^2+1)
    void div(Complex &c, const Complex &a, const Complex &b);

    /*open and reveal*/
    // vs --> vc
    void open(const vector<Share> &vs, vector<gfp> &vc);

    // secrets, clears
    void reveal(const vector<Share> &vs, vector<gfp> &vc);
    void reveal_and_print(const vector<Share> &vs, vector<gfp> &vc);
    void reveal_and_print(const vector<Share> &vs);
    void reveal_and_print(const vector<Complex> &vs);

    /*inputs*/
    void get_inputs(vector<Share> &inputs);
    void get_inputs_dumy(vector<Share> &inputs);

    // the following apis for testing
    void test_add();
    void test_add_plain();
    void test_mul_plain();
    void test_mul();
    void test_sqr();
    void test_div();

    void test_complex_add();
    void test_complex_sub();
    void test_complex_mul_plain();
    void test_complex_mul();
    void test_complex_sqr();
    void test_complex_inv();
    void test_complex_div();
};

#endif