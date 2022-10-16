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

    void mclBnG1_to_gfp(G1_Affine_Coordinates_Plain<T> point, const mclBnG1 &a)
    {
        vector<string> m;
        mclBnG1_to_str_new(m, a);
        cout<<"hhhh"<<m[0].c_str()<<endl;
        cout<<"xxxx"<<m[1].c_str()<<endl;
        cout<<"dddd"<<m[2].c_str()<<endl;
        mpz_class bnx(m[1], 10);
        bigint bn(bnx);
        to_gfp(point.x, bn);
        cout<<"mmmm"<<point.x<<endl;
        mpz_class bnx1(m[2], 10);
        bigint bn1(bnx1);
        to_gfp(point.y, bn1);
        cout<<"uuuuuu"<<point.y<<endl;
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

    void mclBnG2_to_g2_Plain(G2_Affine_Coordinates_Plain<T> &cp, const mclBnG2 &s)
    {

        vector<typename T::clear> vg;
        mclBnG2_to_gfp(vg, s);

        cp.x.real = vg[0];
        cp.x.imag = vg[1];
        cp.y.real = vg[2];
        cp.y.imag = vg[3];
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



    void compute_pk(vector<T> &srk,vector<mclBnFr> &coeff_r, Player &P,
                   typename T::Protocol &protocol, typename T::LivePrep &preprocessing, SubProcessor<T> &processor,
                   typename T::MAC_Check &output)
    {

        
        
        string strp = "52435875175126190479447740508185965837690552500527637822603658699938581184513";
        mpz_class bnx(strp, 10);
        bigint bn(bnx);

        coeff_r.resize(RSIZE);
        vector<mclBnG1> pre_g1(RSIZE);
        mclBnG1 basePoint;
        
        getBasePointG1(basePoint);

        vector <G1_Affine_Coordinates<T> >  r_tmp(255);
        G1Op<T> g1op(P, protocol, preprocessing, processor, output);
        OnlineOp<T> online_op(P, protocol, preprocessing, processor, output);
        for(int i=0;i<255;i++){
            bigint coeff=powerMod(2,i,bn);
            //cout<<"bigbigbig"<<coeff<<endl;
            string coeff_str=to_string(coeff);
            //cout<<"stringinginnig"<<coeff_str.c_str()<<endl;
            str_to_mclBnFr(coeff_r[i],coeff_str);
            //print_mclBnFr(coeff_r[i]);
            mclBnG1_mul(&pre_g1[i], &basePoint, &coeff_r[i]);
            //print_mclBnG1(pre_g1[i]);
            // for(int i=0;i<m.size();i++){
            //     cout<<m[i].c_str() << endl;
            // }
            G1_Affine_Coordinates_Plain<T> point;
            //mclBnG1_to_gfp(point, pre_g1[i]);

            vector<string> m;
            mclBnG1_to_str_new(m, pre_g1[i]);
            // mpz_class bnx(m[1], 10);
            // bigint bn(bnx);
            //point.x=bn;
            online_op.str_to_gfp(point.x, m[1]);
            // mpz_class bnx1(m[2], 10);
            // bigint bn1(bnx1);
            online_op.str_to_gfp(point.y, m[2]);
            //point.y=bn1;

            // cout<<"point.x"<<point.x<< endl;
            // cout<<"point.y"<<point.y<< endl;

            g1op.fixmul_plain_aff(r_tmp[i],srk[i], point);

            // clear m3,m4;
            // online_op.reveal(r_tmp[i].x,m3);
            // cout<<i<<"rtmp"<<m3<<endl;
            // online_op.reveal(r_tmp[i].y,m4);
            // cout<<i<<"rtmp"<<m4<<endl;
            
            // clear srk1;
            // online_op.reveal(srk[i],srk1);

            // cout<<"iiiiiiiiii"<<srk1<<endl;


            T tmp;
            clear c=1;
            online_op.sub_plain_new(tmp, c, srk[i]);
            G1_Affine_Coordinates_Plain<T> o;
            o.x=0;
            o.y=0;
            G1_Affine_Coordinates<T> tmp1;
            g1op.fixmul_plain_aff(tmp1, tmp, o);

            online_op.add_inplace(r_tmp[i].x,tmp1.x);
            online_op.add_inplace(r_tmp[i].y,tmp1.y);


        }

        // cout<<"hhhhhhhhhhhh"<<endl;

        G1_Affine_Coordinates<T> pk;
        bool a=0;

        vector<T> h1(1);
        preprocessing.get_one(DATA_BIT, h1[0]);

        Timer point_add_time;
        point_add_time.start();
        pk=r_tmp[0];
        for(int i=1;i<255;i++){
           
            g1op.add_aff_inplace(pk, r_tmp[i]);

        }

        point_add_time.stop();

        clear pk_x,pk_y;
        online_op.reveal({pk.x},pk_x);
        online_op.reveal({pk.y},pk_y);
        cout<<"public key x" <<pk_x<<endl;
        cout<<"public key y" << pk_y <<endl;



        cout << "point add in G1: " << point_add_time.elapsed() << " seconds" << endl;

#if OK_CODE
        cout << "used triple: " << g2op.UT.UsedTriples << endl;
        cout << "used square: " << g2op.UT.UsedSquares << endl;
        cout << "used bit: " << g2op.UT.UsedBit << endl;
        cout << "used input mask: " << g2op.UT.UsedInputMask << endl;
#endif
    }



void compute_ek(vector<T> &ek,vector<T> &srk,const string &msg, vector<mclBnFr> &coeff_r, Player &P, 
                   typename T::Protocol &protocol, typename T::LivePrep &preprocessing, SubProcessor<T> &processor,
                   typename T::MAC_Check &output)
    {
        
        
        vector<mclBnG2> pre_g2(RSIZE);
        mclBnG2 sig;
        
        mclBnG2_hashAndMapTo(&sig, (const char *)msg.c_str(), msg.size());

        vector <G2_Affine_Coordinates<T> >  sig_tmp(255);
        G2Op<T> g2op(P, protocol, preprocessing, processor, output);
        OnlineOp<T> online_op(P, protocol, preprocessing, processor, output);
        for(int i=0;i<255;i++){
            mclBnG2_mul(&pre_g2[i], &sig, &coeff_r[i]);
            
            G2_Affine_Coordinates_Plain<T> point;

            mclBnG2_to_g2_Plain(point, pre_g2[i]);


            g2op.fixmul_plain_aff(sig_tmp[i],srk[i], point);

            T tmp;
            clear c=1;
            online_op.sub_plain_new(tmp, c, srk[i]);
            G2_Affine_Coordinates_Plain<T> o;
            o.x.real=0;
            o.x.imag=0;
            o.y.real=0;
            o.y.imag=0;
            G2_Affine_Coordinates<T> tmp1;
            g2op.fixmul_plain_aff(tmp1, tmp, o);

            online_op.add_inplace(sig_tmp[i].x.real,tmp1.x.real);
            online_op.add_inplace(sig_tmp[i].x.imag,tmp1.x.imag);
            online_op.add_inplace(sig_tmp[i].y.real,tmp1.y.real);
            online_op.add_inplace(sig_tmp[i].y.imag,tmp1.y.imag);

            // clear m1,m2;
            // online_op.reveal(r_tmp[i].x,m1);
            // cout<<i<<"xxxxxxxxxx"<<m1<<endl;
            // online_op.reveal(r_tmp[i].y,m2);
            // cout<<i<<"yyyyyyyyy"<<m2<<endl;
        }

        // cout<<"hhhhhhhhhhhh"<<endl;

        G2_Affine_Coordinates<T> sigek;
        bool a=0;

        vector<T> h1(1);
        preprocessing.get_one(DATA_BIT, h1[0]);
        
        Timer point_add_time;
        point_add_time.start();
        sigek=sig_tmp[0];
        for(int i=1;i<255;i++){  
            g2op.add_aff_inplace(sigek, sig_tmp[i]);

        }
        ek.resize(2);
        ek[0]=sigek.x.real;
        ek[1]=sigek.x.imag;

        point_add_time.stop();


        cout << "point add time in G2: " << point_add_time.elapsed() << " seconds" << endl;

#if OK_CODE
        cout << "used triple: " << g2op.UT.UsedTriples << endl;
        cout << "used square: " << g2op.UT.UsedSquares << endl;
        cout << "used bit: " << g2op.UT.UsedBit << endl;
        cout << "used input mask: " << g2op.UT.UsedInputMask << endl;
#endif
    }




    void dstb_keygen_new(T &sk_share,Player &P,typename T::Protocol &protocol, typename T::LivePrep &preprocessing, SubProcessor<T> &processor,
                   typename T::MAC_Check &output)
    {
        VSS v(nparty, threshold);
        vector<bls_sk> shs;
        vector<bls_vk> aux;
        v.rnd_secret();
        mclBnFr a=v.get_secret();
        OnlineOp<T> online_op(P, protocol, preprocessing, processor, output);
        vector<T> sk_tmp(P.num_players());
     
        string a1;
        mclBnFr_to_str_new(a1,a);
        //printf("%s",a1.c_str());
        
        mpz_class bnx(a1, 10);
        bigint bn(bnx);

        //clear q=1;

        for(int i=0;i<sk_tmp.size();i++){
            online_op.get_inputs(i,sk_tmp[i], bn);
        }

        sk_share=sk_tmp[0];
        for(int i=1;i<sk_tmp.size();i++){
            online_op.add_inplace(sk_share , sk_tmp[i]);
            
        }
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
