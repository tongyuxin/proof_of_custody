#pragma once

#include "Math/bls12_381.h"
#include "Group.h"

typedef mclBnFr bls_sk;
typedef mclBnG1 bls_vk;
typedef mclBnG2 bls_sigma;

class BLS_
{
    bls_sk sk;

public:
    bls_vk vk;
    bls_sigma sigma;
    unsigned int nparty;
    unsigned int threshold;

    BLS_()
    {
        mclBn_init(MCL_BLS12_381, MCLBN_COMPILED_TIME_VAR);
    }

    BLS_(unsigned int np, unsigned int th) : nparty(np), threshold(th)
    {
        mclBn_init(MCL_BLS12_381, MCLBN_COMPILED_TIME_VAR);
    }

    void set_vk(const bls_vk _vk)
    {
        vk = _vk;
    }

    void set_sigma(const bls_sigma _sigma)
    {
        sigma = _sigma;
    }

    bls_sk get_sk()
    {
        return sk;
    }

    //normal keygen,sign,verify algorithms
    void keygen();
    void sign(const string msg);
    int verify(const bls_sigma _sigma, const string msg);

    //sign with scale*sk on msg
    void sign_scale(const bls_sk scale, const string msg);
    
    void sign_scale_test(mclBnG2 &wu,const bls_sk scale, const string msg, Player &P);

    //distributed keygen,sign algorithms
    void dstb_keygen(Player &P);
};

//
//
//
//
//

#include "vss.h"
#include "Tools/Commit.h"
#include "Math/Lagrange.h"

template <class T>
class BLS : public BLS_
{
    typedef typename T::clear clear;

public:
    using BLS_::BLS_;

    void mclBnG2_to_gfp(vector<typename T::clear> &out, const mclBnG2 &a)
    {
        vector<string> str;
        mclBnG2_to_str(str, a);
        out.resize(str.size());

        for (int i = 0; i < out.size(); i++)
        {
            mpz_class bnx(str[i], 10);
            bigint bn(bnx);
            to_gfp(out[i], bn);
        }
    }

    void mclBnG2_to_Complex_Plain(vector<Complex_Plain<T>> &cp, const bls_sigma &s)
    {
        cp.resize(2);

        vector<typename T::clear> vg;
        mclBnG2_to_gfp(vg, s);

        cp[0].real = vg[0];
        cp[0].imag = vg[1];
        cp[1].real = vg[2];
        cp[1].imag = vg[3];
    }

    void dstb_sign(G2_Affine_Coordinates<T> &out, const string msg, Player &P,
                   typename T::Protocol &protocol, typename T::LivePrep &preprocessing, SubProcessor<T> &processor,
                   typename T::MAC_Check &output)
    {
        Timer point_add_time;

        mclBnFr lag_coeff;
        get_lagrange_coeff(lag_coeff, P.num_players(), P.my_num() + 1);

        sign_scale(lag_coeff, msg);

        vector<G2_Affine_Coordinates<T>> ac(P.num_players());
        vector<vector<Complex_Plain<T>>> s(P.num_players(), vector<Complex_Plain<T>>(2));

        mclBnG2_to_Complex_Plain(s[P.my_num()], sigma);

        G2Op<T> g2op(P, protocol, preprocessing, processor, output);
        OnlineOp<T> online_op(P, protocol, preprocessing, processor, output);

        // cout<<"getinput before"<<endl;
        // P.comm_stats.print();

        for (int i = 0; i < ac.size(); i++)  
        {

            online_op.get_inputs(i, ac[i].x.real, s[i][0].real);
            online_op.get_inputs(i, ac[i].x.imag, s[i][0].imag);
            online_op.get_inputs(i, ac[i].y.real, s[i][1].real);
            online_op.get_inputs(i, ac[i].y.imag, s[i][1].imag);

            // cout<<"getinput one"<<endl;
            // P.comm_stats.print();
            //g2op.get_inputs(i, ac[i].x, s[i][0]);
           // g2op.get_inputs(i, ac[i].y, s[i][1]);
        }

        cout<<"aff add before"<<endl;
        P.comm_stats.print();

        point_add_time.start();
        out = ac[0];

        for (int i = 1; i < ac.size(); i++) 
        {
            g2op.add_aff_inplace(out, ac[i]);
        }

        cout<<"aff add after"<<endl;
        P.comm_stats.print();
        point_add_time.stop();

        cout << "point addition time (including part of the offline time): " << point_add_time.elapsed() << " seconds" << endl;

#if OK_CODE
        cout << "used triple: " << g2op.UT.UsedTriples << endl;
        cout << "used square: " << g2op.UT.UsedSquares << endl;
        cout << "used bit: " << g2op.UT.UsedBit << endl;
        cout << "used input mask: " << g2op.UT.UsedInputMask << endl;
#endif
    }

    void dstb_sign_test(vector<bigint> &out, const string msg, Player &P,
                   typename T::Protocol &protocol, typename T::LivePrep &preprocessing, SubProcessor<T> &processor,
                   typename T::MAC_Check &output)
    {
        Timer point_add_time;


        mclBnFr lag_coeff;
        get_lagrange_coeff(lag_coeff, P.num_players(), P.my_num() + 1);
        
        mclBnG2 sigma_plain;
        sign_scale_test(sigma_plain,lag_coeff, msg, P);

        cout<<"signature plain: "<<endl;
        print_mclBnG2(sigma_plain);

        vector<string> sigma_str;
        mclBnG2_to_str(sigma_str,sigma_plain);
        out.resize(4);
        for(int i = 0; i < sigma_str.size(); i++)
        {
            mpz_class bnx(sigma_str[i], 10);
            bigint bn(bnx);
            //to_gfp(out[i], bn);
            out[i]=bn;
        }

        // out.x.real = ac[0].real;
        // out.x.imag = ac[0].imag;
        // out.y.real = ac[1].real;
        // out.y.imag = ac[1].imag;



#if OK_CODE
        cout << "used triple: " << g2op.UT.UsedTriples << endl;
        cout << "used square: " << g2op.UT.UsedSquares << endl;
        cout << "used bit: " << g2op.UT.UsedBit << endl;
        cout << "used input mask: " << g2op.UT.UsedInputMask << endl;
#endif
    }



};
