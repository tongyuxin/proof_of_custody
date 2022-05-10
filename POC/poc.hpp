#include "poc.h"

template <class T>
void POC<T>::poc_setup(BLS<T> &bls, Player &P)
{
    bls.dstb_keygen(P);
}

template <class T>
void POC<T>::poc_compute_ephem_key(G2_Affine_Coordinates<T> &out, BLS<T> &bls, const string &msg, int online_num, Player &P, Config_Info &CI)
{
    bls.dstb_sign(out, msg, P, protocol, preprocessing, processor, output);
}

// pre_key = (s_1, s_0^2, s_1^3, s_0^4,...,)
template <class T>
void POC<T>::poc_compute_custody_bit_offline(vector<T> &pre_key, const vector<T> &keys, int online_num, Player &P, Config_Info &CI)
{
    if (keys.size() != 2)
    {
        throw invalid_length();
    }

    OnlineOp<T> online_op(P, protocol, preprocessing, processor, output);

    pre_key.resize(CHUNK_NUM);
    T s0_square, s1_square;

    online_op.sqr(s0_square, keys[0]);
    online_op.sqr(s1_square, keys[1]);

    pre_key[0] = keys[1];
    pre_key[1] = s0_square;

    for (int i = 2; i < CHUNK_NUM; i++)
    {
        if (i % 2 == 0)
        {
            online_op.mul(pre_key[i], pre_key[i - 2], s1_square);
        }
        else
        {
            online_op.mul(pre_key[i], pre_key[i - 2], s0_square);
        }
    }

    //*RC*// cout << "used triple: " << online_op.UT.UsedTriples << endl;
    //*RC*// cout << "used square: " << online_op.UT.UsedSquares << endl;
    //*RC*// cout << "used bit: " << online_op.UT.UsedBit << endl;
    //*RC*// cout << "used input mask: " << online_op.UT.UsedInputMask << endl;
}

template <class T>
int POC<T>::poc_compute_custody_bit_online(const vector<T> pre_key, const vector<clear> &msg, int online_num, Player &P, Config_Info &CI)
{
    if (msg.size() != CHUNK_NUM)
    {
        throw invalid_length();
    }
    OnlineOp<T> online_op(P, protocol, preprocessing, processor, output);

    T uhf_out;

    online_op.mul_plain(uhf_out, pre_key[0], msg[1]);
    online_op.add_plain_inplace(uhf_out, msg[0]);

    T tmp;

    for (int i = 2; i < msg.size(); i++)
    {
        online_op.mul_plain(tmp, pre_key[i - 1], msg[i]);
        online_op.add_inplace(uhf_out, tmp);
    }

    online_op.add_inplace(uhf_out, pre_key.back());

    int res;
    res = online_op.legendre_prf(pre_key[0], uhf_out);

    //*RC*// cout << "used triple: " << online_op.UT.UsedTriples << endl;
    //*RC*// cout << "used square: " << online_op.UT.UsedSquares << endl;
    //*RC*// cout << "used bit: " << online_op.UT.UsedBit << endl;
    //*RC*// cout << "used input mask: " << online_op.UT.UsedInputMask << endl;

    return res;
}

template <class T>
int POC<T>::poc_compute_custody_bit(const vector<T> &keys, const vector<clear> &msg, int online_num, Player &P, Config_Info &CI)
{
    if (keys.size() != 2)
    {
        throw invalid_length();
    }

    OnlineOp<T> online_op(P, protocol, preprocessing, processor, output);

    T uhf_out;
    int res;
    online_op.uhf(uhf_out, keys[0], msg, CHUNK_NUM);
    res = online_op.legendre_prf(keys[1], uhf_out);

    //*RC*// cout << "used triple: " << online_op.UT.UsedTriples << endl;
    //*RC*// cout << "used square: " << online_op.UT.UsedSquares << endl;
    //*RC*// cout << "used bit: " << online_op.UT.UsedBit << endl;
    //*RC*// cout << "used input mask: " << online_op.UT.UsedInputMask << endl;

    return res;
}

template <class T>
void POC<T>::shared_rand_bits_phase_one_new(vector<T> &shared_bits,vector<bigint> &local_bits,vector<bigint> &sigma_bits,vector<bigint> &x_bits_plain,  int online_num, Player &P,
                                        Config_Info &CI)
{
    Timer daBits_phase_one;
    daBits_phase_one.start();
    local_bits.resize(2 * PSIZE + SEC);
    shared_bits.resize(2 * PSIZE + SEC);

    vector<uint8_t> rnd(101); // (2*381+40)/8
    //*RC*// P.G.get_random_bytes(rnd);
    PRNG prng;
    prng.ReSeed();
    prng.get_octets((octet *)rnd.data(), sizeof(uint8_t) * rnd.size());

    for (int i = 0; i < local_bits.size(); i++)
    {
        local_bits[i] = (rnd[i / 8] >> (i % 8)) & 1;
    }

    OnlineOp<T> online_op(P, protocol, preprocessing, processor, output);
    vector<T> sbit(P.num_players());
    for (int j = 0; j < 2 * PSIZE + SEC; j++)
    {
        for (int i = 0; i < sbit.size(); i++)
        {
            clear tmp;
            tmp = local_bits[j];
            //sbit[i].set_player(P.whoami());
            //if (i == P.whoami()) {
            online_op.get_inputs(i, sbit[i], tmp);                                
            //}
        }
        online_op.KXOR(shared_bits[j], sbit, sbit.size()); 
        //shared_bits_tmp[j]=shared_bits[j];                        
    }

    //consistency check
    //r
    vector<T> r_bits(SEC+2);
    vector<T> tmp1(1);
    for (int i = 0; i < SEC+2; i++)
    {
      preprocessing.get_one(DATA_BIT, tmp1[0]);
      // getTuples(tmp, BIT);
      r_bits[i] = tmp1[0];
    }

    // vector<uint8_t> rnd1(96); // 2*381/8
    // //*RC*// P.G.get_random_bytes(rnd);
    // PRNG prng1;
    // prng1.ReSeed();
    // prng1.get_octets((octet *)rnd1.data(), sizeof(uint8_t) * rnd1.size());
    //sigma 
    
    sigma_bits.resize(2 * PSIZE);
    vector<T> sigma_b(2 * PSIZE);
    vector<T> tmp_sigma(1);

    for (int i = 0; i < sigma_b.size(); i++)
    {
        preprocessing.get_one(DATA_BIT, tmp_sigma[0]);
        sigma_b[i]=tmp_sigma[0];
        //sigma_bits[i] =(rnd1[i / 8] >> (i % 8)) & 1  ;
    }
    vector<clear> sigma_reveal;
    online_op.reveal({sigma_b},sigma_reveal);
    for (int i = 0; i < sigma_bits.size(); i++)
    {
        to_bigint(sigma_bits[i], sigma_reveal[i]);
    }
    

    //x
    vector<T> x_bits(SEC);
    vector<T> tmp3(SEC);
    vector<T> tmp2(2 * PSIZE);
    bigint TWO(2);
    
    for(int i=0;i<SEC;i++)
    {
        for(int j=0;j< 2 * PSIZE;j++)
        {
            online_op.mul_plain(tmp2[j],shared_bits[j],sigma_bits[j]);
            online_op.add_inplace(x_bits[i],tmp2[j]);
        }
        online_op.add_inplace(x_bits[i],shared_bits[2 * PSIZE+i]);
        
       
        vector<T> tmp4(SEC+2);
        
        for(int j=0;j<SEC+2;j++)
        {
            bigint m= (TWO << j);
            online_op.mul_plain(tmp4[j],r_bits[j], m);
            online_op.add_inplace(tmp3[i],tmp4[j]);
        }
        online_op.mul_plain(tmp3[i],tmp3[i],2);
        online_op.add_inplace(x_bits[i],tmp3[i]);
    }

    x_bits_plain.resize(SEC);
    vector<bigint> x_bits_tmp(SEC);
    vector<clear> x_tmp;
    online_op.reveal({x_bits},x_tmp);
    for(int i=0;i<SEC;i++){
        to_bigint(x_bits_tmp[i], x_tmp[i]);
        x_bits_plain[i]=x_bits_tmp[i] % 2;
    }

    daBits_phase_one.stop();
    cout << "Time required for the first phase of daBits generation and verification: " << daBits_phase_one.elapsed() << " seconds" << endl;

    //*RC*// cout << "used triple: " << online_op.UT.UsedTriples << endl;
    //*RC*// cout << "used square: " << online_op.UT.UsedSquares << endl;
    //*RC*// cout << "used bit: " << online_op.UT.UsedBit << endl;
    //*RC*// cout << "used input mask: " << online_op.UT.UsedInputMask << endl;
}



template <class T>
void POC<T>::shared_rand_bits_phase_one(vector<T> &shared_bits, vector<bigint> &local_bits, int online_num, Player &P,
                                        Config_Info &CI)
{
    Timer phase_one;
    phase_one.start();

    local_bits.resize(2 * PSIZE);
    shared_bits.resize(2 * PSIZE);

    vector<uint8_t> rnd(96); // 2*381/8
    //*RC*// P.G.get_random_bytes(rnd);
    PRNG prng;
    prng.ReSeed();
    prng.get_octets((octet *)rnd.data(), sizeof(uint8_t) * rnd.size());

    for (int i = 0; i < local_bits.size(); i++)
    {
        local_bits[i] = (rnd[i / 8] >> (i % 8)) & 1;
    }

    OnlineOp<T> online_op(P, protocol, preprocessing, processor, output);

    vector<T> sbit(P.num_players());
    for (int j = 0; j < 2 * PSIZE; j++)
    {
        for (int i = 0; i < sbit.size(); i++)
        {
            clear tmp;
            // tmp.assign(local_bits[j]);
            tmp = local_bits[j];
            //sbit[i].set_player(P.whoami());
            //if (i == P.whoami()) {
            online_op.get_inputs(i, sbit[i], tmp);                                
            //}
        }
        online_op.KXOR(shared_bits[j], sbit, sbit.size());                              
    }

    //*RC*// cout << "used triple: " << online_op.UT.UsedTriples << endl;
    //*RC*// cout << "used square: " << online_op.UT.UsedSquares << endl;
    //*RC*// cout << "used bit: " << online_op.UT.UsedBit << endl;
    //*RC*// cout << "used input mask: " << online_op.UT.UsedInputMask << endl;
    phase_one.stop();
    cout << "Time required for the first phase of random shared bit generation: " << phase_one.elapsed() << " seconds" << endl;
}
template <class T>
void POC<T>::decompose_and_reveal(
    vector<bigint> &reveal_bits, const vector<T> &keys, const vector<T> &shared_bits, 
    int online_num, Player &P, Config_Info &CI)
{
    if (keys.size() != 2 || shared_bits.size() != 2 * PSIZE)
    {
        throw invalid_length();
    }

    PRINT_DEBUG_INFO();
    reveal_bits.resize(2 * PSIZE);

    OnlineOp<T> online_op(P, protocol, preprocessing, processor, output);
    PRINT_DEBUG_INFO();

    vector<T> deco_bits, tmp;
    PRINT_DEBUG_INFO();
    online_op.A2B(deco_bits, keys[0]);
    PRINT_DEBUG_INFO();
    online_op.A2B(tmp, keys[1]);
    PRINT_DEBUG_INFO();
    deco_bits.insert(deco_bits.end(), tmp.begin(), tmp.end());

    PRINT_DEBUG_INFO();
    for (int i = 0; i < shared_bits.size(); i++)
    {
        online_op.XOR_inplace(deco_bits[i], shared_bits[i]);
    }

    PRINT_DEBUG_INFO();
    vector<clear> out;
    online_op.reveal(deco_bits, out);
    for (int i = 0; i < reveal_bits.size(); i++)
    {
        to_bigint(reveal_bits[i], out[i]);
    }
    PRINT_DEBUG_INFO();

    //*RC*// cout << "used triple: " << online_op.UT.UsedTriples << endl;
    //*RC*// cout << "used square: " << online_op.UT.UsedSquares << endl;
    //*RC*// cout << "used bit: " << online_op.UT.UsedBit << endl;
    //*RC*// cout << "used input mask: " << online_op.UT.UsedInputMask << endl;
}



template <class T>
void POC<T>::decompose_and_reveal_new(
    vector<bigint> &reveal_bits, const vector<T> &keys, const vector<T> &shared_bits, 
    int online_num, Player &P, Config_Info &CI)
{
    if (keys.size() != 2 || shared_bits.size() != 2 * PSIZE+SEC)
    {
        throw invalid_length();
    }

    PRINT_DEBUG_INFO();
    reveal_bits.resize(2 * PSIZE);

    OnlineOp<T> online_op(P, protocol, preprocessing, processor, output);
    PRINT_DEBUG_INFO();

    vector<T> deco_bits, tmp;
    PRINT_DEBUG_INFO();
    online_op.A2B(deco_bits, keys[0]);
    PRINT_DEBUG_INFO();
    online_op.A2B(tmp, keys[1]);
    PRINT_DEBUG_INFO();
    deco_bits.insert(deco_bits.end(), tmp.begin(), tmp.end());

    PRINT_DEBUG_INFO();
    for (int i = 0; i < 2 * PSIZE; i++)
    {
        online_op.XOR_inplace(deco_bits[i], shared_bits[i]);
    }

    PRINT_DEBUG_INFO();
    vector<clear> out;
    online_op.reveal(deco_bits, out);
    for (int i = 0; i < reveal_bits.size(); i++)
    {
        to_bigint(reveal_bits[i], out[i]);
    }
    PRINT_DEBUG_INFO();

    //*RC*// cout << "used triple: " << online_op.UT.UsedTriples << endl;
    //*RC*// cout << "used square: " << online_op.UT.UsedSquares << endl;
    //*RC*// cout << "used bit: " << online_op.UT.UsedBit << endl;
    //*RC*// cout << "used input mask: " << online_op.UT.UsedInputMask << endl;
}

template <class T>
void POC<T>::shared_rand_bits_phase_two_new(
    vector<T> &shared_bits, const vector<bigint> &local_bits,const vector<bigint> &sigma_bits,
    const vector<bigint> &x_bits, int online_num, Player &P,Config_Info &CI)
{
    Timer daBits_phase_two;
    daBits_phase_two.start();
    shared_bits.resize(2 * PSIZE+SEC);
    OnlineOp<T> online_op(P, protocol, preprocessing, processor, output);

    vector<T> sbit(P.num_players());

    for (int j = 0; j < 2 * PSIZE + SEC; j++)
    {

        for (int i = 0; i < sbit.size(); i++)
        {
            clear tmp;
            // tmp.assign(local_bits[j]);
            tmp = local_bits[j];
            //sbit[i].set_player(P.whoami());
            //if (i == P.whoami()) {
            online_op.get_inputs(i, sbit[i], tmp);
            //}
        }
        online_op.KXOR(shared_bits[j], sbit, sbit.size());
    }


    //consistency Check

    //t
    vector<T> t_bits(SEC+2);
    vector<T> tmp1(1);
    for (int i = 0; i < SEC+2; i++)
    {
      preprocessing.get_one(DATA_BIT, tmp1[0]);
      // getTuples(tmp, BIT);
      t_bits[i] = tmp1[0];
    }
    //y
    vector<T> y_bits(SEC);
    vector<T> tmp2(2 * PSIZE);
    vector<T> tmp3(SEC);
    bigint TWO(2);
    for(int i=0;i<SEC;i++)
    {
        for(int j=0;j< 2 * PSIZE;j++){
            online_op.mul_plain(tmp2[j],shared_bits[j],sigma_bits[j]);
            online_op.add_inplace(y_bits[i],tmp2[j]);
        }
        online_op.add_inplace(y_bits[i],shared_bits[2 * PSIZE+i]);
        
       
        vector<T> tmp4(SEC+2);
        for(int j=0;j<SEC+2;j++)
        {
            bigint m = TWO << j;
            online_op.mul_plain(tmp4[j],t_bits[j],m );
            online_op.add_inplace(tmp3[i],tmp4[j]);
        }
        online_op.mul_plain(tmp3[i],tmp3[i],2);
        online_op.add_inplace(y_bits[i],tmp3[i]);
    }

    vector<bigint> y_bits_plain(SEC);
    cout<<"Perform consistency check of daBits"<<endl;
    
    vector<clear> y_tmp;
    online_op.reveal({y_bits},y_tmp);

    for(int i=0;i<SEC;i++)
    {
        to_bigint(y_bits_plain[i], y_tmp[i]);
        if(x_bits[i]  != y_bits_plain[i] % 2)
        {
            cout<<x_bits[i]<<endl;
            cout<<y_bits_plain[i] % 2<<endl;
            throw bad_dabits_value();

        }
    }
    cout<<"The consistency check passed!!!"<<endl;
    daBits_phase_two.stop();
    cout << "Time required for the second phase of daBits generation and verification: " << daBits_phase_two.elapsed() << " seconds" << endl;

    //*RC*// cout << "used triple: " << online_op.UT.UsedTriples << endl;
    //*RC*// cout << "used square: " << online_op.UT.UsedSquares << endl;
    //*RC*// cout << "used bit: " << online_op.UT.UsedBit << endl;
    //*RC*// cout << "used input mask: " << online_op.UT.UsedInputMask << endl;
}

template <class T>
void POC<T>::shared_rand_bits_phase_two(
    vector<T> &shared_bits, const vector<bigint> &local_bits, int online_num, Player &P,
    Config_Info &CI)
{
    Timer phase_two;
    phase_two.start();
    shared_bits.resize(2 * PSIZE);
    OnlineOp<T> online_op(P, protocol, preprocessing, processor, output);

    vector<T> sbit(P.num_players());
    for (int j = 0; j < 2 * PSIZE; j++)
    {
        for (int i = 0; i < sbit.size(); i++)
        {
            clear tmp;
            // tmp.assign(local_bits[j]);
            tmp = local_bits[j];
            //sbit[i].set_player(P.whoami());
            //if (i == P.whoami()) {
            online_op.get_inputs(i, sbit[i], tmp);
            //}
        }
        online_op.KXOR(shared_bits[j], sbit, sbit.size());
    }

    phase_two.stop();
    cout << "Time required for the second phase of random shared bit generation: " << phase_two.elapsed() << " seconds" << endl;


    //*RC*// cout << "used triple: " << online_op.UT.UsedTriples << endl;
    //*RC*// cout << "used square: " << online_op.UT.UsedSquares << endl;
    //*RC*// cout << "used bit: " << online_op.UT.UsedBit << endl;
    //*RC*// cout << "used input mask: " << online_op.UT.UsedInputMask << endl;
}

template <class T>
void POC<T>::xor_and_combine(
    vector<T> &keys, const vector<T> &shared_bits, const vector<bigint> &reveal_bits,
    int online_num, Player &P, Config_Info &CI)
{
    if (shared_bits.size() != 3 * QSIZE || reveal_bits.size() != 3 * QSIZE)
    {
        throw invalid_length();
    }

    keys.resize(3);
    OnlineOp<T> online_op(P, protocol, preprocessing, processor, output);

    vector<T> tmp(shared_bits.size());
    for (int i = 0; i < shared_bits.size(); i++)
    {
        online_op.XOR_plain(tmp[i], shared_bits[i], reveal_bits[i]);
    }

    vector<T> out;
    for (int i = 0; i < keys.size(); i++)
    {
        out.insert(out.begin(), tmp.begin() + i * QSIZE, tmp.begin() + (i + 1) * QSIZE);
        online_op.B2A(keys[i], out, QSIZE);
        out.clear();
    }
    tmp.clear();

    //*RC*// cout << "used triple: " << online_op.UT.UsedTriples << endl;
    //*RC*// cout << "used square: " << online_op.UT.UsedSquares << endl;
    //*RC*// cout << "used bit: " << online_op.UT.UsedBit << endl;
    //*RC*// cout << "used input mask: " << online_op.UT.UsedInputMask << endl;
}



template <class T>
void POC<T>::xor_and_combine_new(
    vector<T> &keys, const vector<T> &shared_bits, const vector<bigint> &reveal_bits,
    int online_num, Player &P, Config_Info &CI)
{
    if (shared_bits.size() != 3 * QSIZE + SEC || reveal_bits.size() != 3 * QSIZE)
    {
        throw invalid_length();
    }

    keys.resize(3);
    OnlineOp<T> online_op(P, protocol, preprocessing, processor, output);

    vector<T> tmp(reveal_bits.size());
    for (int i = 0; i < reveal_bits.size(); i++)
    {
        online_op.XOR_plain(tmp[i], shared_bits[i], reveal_bits[i]);
    }

    vector<T> out;
    for (int i = 0; i < keys.size(); i++)
    {
        out.insert(out.begin(), tmp.begin() + i * QSIZE, tmp.begin() + (i + 1) * QSIZE);
        online_op.B2A(keys[i], out, QSIZE);
        out.clear();
    }
    tmp.clear();

    //*RC*// cout << "used triple: " << online_op.UT.UsedTriples << endl;
    //*RC*// cout << "used square: " << online_op.UT.UsedSquares << endl;
    //*RC*// cout << "used bit: " << online_op.UT.UsedBit << endl;
    //*RC*// cout << "used input mask: " << online_op.UT.UsedInputMask << endl;
}

template <class T>
void POC<T>::poc_compute_custody_bit_offline_2primes(
    vector<T> &pre_key, const vector<T> &keys, int online_num, Player &P, Config_Info &CI)
{
    if (keys.size() != 3)
    {
        throw bad_value();
    }

    OnlineOp<T> online_op(P, protocol, preprocessing, processor, output);
    PRINT_DEBUG_INFO();

    pre_key.resize(CHUNK_NUM);
    T s0_cube, s1_cube, s2_cube;
    online_op.sqr(s0_cube, keys[0]);
    online_op.mul_inplace(s0_cube, keys[0]);

    online_op.sqr(s1_cube, keys[1]);
    online_op.mul_inplace(s1_cube, keys[1]);

    online_op.sqr(s2_cube, keys[2]);
    pre_key[1] = s2_cube;
    online_op.mul_inplace(s2_cube, keys[2]);

    pre_key[0] = keys[1];
    pre_key[2] = s0_cube;
    // s1, s2^2, s0^3, s1^4, s2^5, s0^6, constant not considered

    PRINT_DEBUG_INFO();
    for (int i = 3; i < CHUNK_NUM; i++)
    {
        if (i % 1000 == 0)
            PRINT_DEBUG_INFO();
        if (i % 3 == 0)
        {
            online_op.mul(pre_key[i], pre_key[i - 3], s1_cube);
        }
        else if (i % 3 == 1)
        {
            online_op.mul(pre_key[i], pre_key[i - 3], s2_cube);
        }
        else if (i % 3 == 2)
        {
            online_op.mul(pre_key[i], pre_key[i - 3], s0_cube);
        }
    }

    //*RC*// cout << "used triple: " << online_op.UT.UsedTriples << endl;
    //*RC*// cout << "used square: " << online_op.UT.UsedSquares << endl;
    //*RC*// cout << "used bit: " << online_op.UT.UsedBit << endl;
    //*RC*// cout << "used input mask: " << online_op.UT.UsedInputMask << endl;
}

template <class T>
int POC<T>::poc_compute_custody_bit_online_2primes(
    const vector<T> pre_key, const T key, const vector<clear> &msg, int online_num, Player &P,
    Config_Info &CI)
{
    if (msg.size() != CHUNK_NUM)
    {
        throw bad_value();
    }
    OnlineOp<T> online_op(P, protocol, preprocessing, processor, output);

    T uhf_out;

    online_op.mul_plain(uhf_out, pre_key[0], msg[1]);
    online_op.add_plain_inplace(uhf_out, msg[0]); // m[0] + m[1]*s1

    T tmp;

    for (int i = 2; i < msg.size(); i++)
    {
        online_op.mul_plain(tmp, pre_key[i - 1], msg[i]);
        online_op.add_inplace(uhf_out, tmp);
    }

    online_op.add_inplace(uhf_out, pre_key.back());

    int res = 0;
    for (int i = 0; i < 10; i++)
    {
        clear count(i);
        T in;
        online_op.add_plain(in, uhf_out, count);
        res &= online_op.legendre_prf(key, in);
    }

    //*RC*// cout << "used triple: " << online_op.UT.UsedTriples << endl;
    //*RC*// cout << "used square: " << online_op.UT.UsedSquares << endl;
    //*RC*// cout << "used bit: " << online_op.UT.UsedBit << endl;
    //*RC*// cout << "used input mask: " << online_op.UT.UsedInputMask << endl;

    return res;
}
