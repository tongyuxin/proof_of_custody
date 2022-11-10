#include "poc.h"

template <class T>
void POC<T>::poc_setup(BLS<T> &bls, Player &P)
{
    bls.dstb_keygen(P);
}

template <class T>
void POC<T>::poc_setup_new(T &sk_share, BLS<T> &bls, Player &P)
{
    bls.dstb_keygen_new(sk_share, P, protocol, preprocessing, processor, output);
}

template <class T>
void POC<T>::poc_compute_ephem_key(G2_Affine_Coordinates<T> &out, BLS<T> &bls, const string &msg, int online_num, Player &P, Config_Info &CI)
{
    bls.dstb_sign(out, msg, P, protocol, preprocessing, processor, output);
}

template <class T>
void POC<T>::poc_compute_public_key(vector<T> &srk, vector<mclBnFr> &coeff_r, BLS<T> &bls, int online_num, Player &P, Config_Info &CI)
{
    bls.compute_pk(srk, coeff_r, P, protocol, preprocessing, processor, output);
}

template <class T>
void POC<T>::poc_compute_ek(vector<T> &ek,vector<T> &srk, BLS<T> &bls, int online_num, Player &P, Config_Info &CI,const string &msg,vector<mclBnFr> &coeff_r)
{
    bls.compute_ek(ek,srk, msg, coeff_r, P, protocol, preprocessing, processor, output);
}

template <class T>
void POC<T>::poc_compute_ephem_key_test(vector<bigint> &out, BLS<T> &bls, const string &msg, int online_num, Player &P, Config_Info &CI)
{
    //bls.dstb_sign_test1(out, msg, P, protocol, preprocessing, processor, output);
    bls.dstb_sign_test(out, msg, P, protocol, preprocessing, processor, output);
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
void POC<T>::shared_rand_bits_for_pk_phase_one_new(vector<T> &shared_bits,vector<bigint> &local_bits,vector<bigint> &sigma_bits,vector<bigint> &x_bits_plain,  int online_num, Player &P,
                                        Config_Info &CI)
{
    Timer daBits_phase_one;
    daBits_phase_one.start();
    local_bits.resize(RSIZE + SEC);
    shared_bits.resize(RSIZE + SEC);

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

    // cout<<"before randbit one"<<endl;
    // P.comm_stats.print();
    for (int j = 0; j < RSIZE + SEC; j++)
    {
        for (int i = 0; i < sbit.size(); i++)
        {
            clear tmp;
            tmp = local_bits[j];
            //sbit[i].set_player(P.whoami());
            //if (i == P.whoami()) {
            online_op.get_inputs(i, sbit[i], tmp);
            
            // cout<<"after get input"<<endl;
            // P.comm_stats.print();                          
            //}
        }
        online_op.KXOR(shared_bits[j], sbit, sbit.size()); 

        // cout<<"after xor"<<endl;
        // P.comm_stats.print();  
        //shared_bits_tmp[j]=shared_bits[j];                        
    }

    // cout<<"after randbit one"<<endl;
    // P.comm_stats.print();

    //consistency check
    //r
    vector<T> r_bits(SEC+2);
    vector<T> tmp1(1);

    // cout<<"before get one"<<endl;
    // P.comm_stats.print();
    for (int i = 0; i < SEC+2; i++)
    {
      preprocessing.get_one(DATA_BIT, tmp1[0]);
      // getTuples(tmp, BIT);
      r_bits[i] = tmp1[0];
    }

    // cout<<"after get one"<<endl;
    // P.comm_stats.print();


    // vector<uint8_t> rnd1(96); // 2*381/8
    // //*RC*// P.G.get_random_bytes(rnd);
    // PRNG prng1;
    // prng1.ReSeed();
    // prng1.get_octets((octet *)rnd1.data(), sizeof(uint8_t) * rnd1.size());
    
    
    
    
    //sigma 
    
    // sigma_bits.resize(2 * PSIZE);
    // vector<T> sigma_b(2 * PSIZE);
    // vector<T> tmp_sigma(1);

    // for (int i = 0; i < sigma_b.size(); i++)
    // {
    //     preprocessing.get_one(DATA_BIT, tmp_sigma[0]);
    //     sigma_b[i]=tmp_sigma[0];
    //     //sigma_bits[i] =(rnd1[i / 8] >> (i % 8)) & 1  ;
    // }
    // vector<clear> sigma_reveal;
    // online_op.reveal({sigma_b},sigma_reveal);
    // for (int i = 0; i < sigma_bits.size(); i++)
    // {
    //     to_bigint(sigma_bits[i], sigma_reveal[i]);
    // }

    //sigma 
    vector<bigint> tmp_sigma(RSIZE);
    vector<uint8_t> rnd1(96); // 2*381/8
    //*RC*// P.G.get_random_bytes(rnd);
    PRNG prng1;
    prng1.ReSeed();
    prng1.get_octets((octet *)rnd1.data(), sizeof(uint8_t) * rnd1.size());
    for (int i = 0; i < tmp_sigma.size(); i++)
    {

        tmp_sigma[i] =(rnd1[i / 8] >> (i % 8)) & 1  ;
    }
    stringstream ss;
    for(int i=0; i < tmp_sigma.size(); i++){
        ss << tmp_sigma[i];
    }
    string s(ss.str());
    // for(int i=0;i<s.length();i++){
    //     cout<<s[i];
    // }

    octetStream os_comm;
    octetStream os_open;
    octetStream os_data(s.size(), (const octet *)s.data());
    Commit(os_comm, os_open, os_data, 1); 
    // send my commitment

    // cout<<"before send all"<<endl;
    // P.comm_stats.print();
    P.send_all(os_comm); 
    // cout<<"after send all"<<endl;
    // P.comm_stats.print();

    // receive other commitments
    vector<octetStream> os_comms(P.num_players()); // CommAux(P.num_players());
    P.receive_all(os_comms); 

    // cout<<"after receive all"<<endl;
    // P.comm_stats.print();

    // send my opening
    P.send_all(os_open);

    // cout<<"after send all open"<<endl;
    // P.comm_stats.print();

    // receive other opening
    vector<octetStream> os_opens(P.num_players()); // OpenAux(P.num_players());
    P.receive_all(os_opens);

    // cout<<"after receive all open"<<endl;
    // P.comm_stats.print();


    
    octetStream os_ss;
    vector<string> sigma_share(P.num_players());
    sigma_share[P.my_num()]= s;
    for(int i=0; i<P.num_players();i++){
        if(i!=P.my_num()){
            bool res=Open(os_ss, os_comms[i],os_opens[i],1);
            if(!res){
                throw bad_dabits_value();
            }
            sigma_share[i].assign((char *)os_ss.get_data(), os_ss.get_length());
        }
    }


    sigma_bits.resize(RSIZE);


    // vector<clear> neww;
    // neww.resize(762);
    vector<vector<bigint> > sigma_share_tmp(P.num_players(),vector<bigint> (RSIZE));
    for(int i=0; i<RSIZE;i++){
        for(int j=0;j<P.num_players();j++){
            string s_tmp;
            s_tmp.resize(1);
            s_tmp[0]=sigma_share[j][i];
            mpz_class bnx1(s_tmp, 2);
            bigint bn1(bnx1);
            sigma_share_tmp[j][i]=bn1;
            sigma_bits[i]=sigma_bits[i] ^ sigma_share_tmp[j][i];
        }

        // cout<< sigma_bits[i];
    }

    

    //x
    vector<T> x_bits(SEC);
    vector<T> tmp3(SEC);
    vector<T> tmp2(RSIZE);
    bigint TWO(2);
    
    for(int i=0;i<SEC;i++)
    {
        for(int j=0;j< RSIZE;j++)
        {
            online_op.mul_plain(tmp2[j],shared_bits[j],sigma_bits[j]);
            online_op.add_inplace(x_bits[i],tmp2[j]);
        }
        online_op.add_inplace(x_bits[i],shared_bits[RSIZE+i]);
        
       
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
    
    // cout<<"before reveal"<<endl;
    // P.comm_stats.print();
    online_op.reveal({x_bits},x_tmp);
    
    // cout<<"after reveal"<<endl;
    // P.comm_stats.print();
    for(int i=0;i<SEC;i++){
        to_bigint(x_bits_tmp[i], x_tmp[i]);
        x_bits_plain[i]=x_bits_tmp[i] % 2;
    }

     cout<<"Perform bit check of daBits"<<endl;
    
    for(int i=0;i<RSIZE;i++)
    {
        clear v=1;
        T u;
        online_op.sub_plain_new(u,v,shared_bits[i]);
        online_op.mul_inplace(u,shared_bits[i]);
        clear c;
        online_op.reveal(u,c);
        if(c != 0)
        {
            throw bad_dabits_value();
        }
    }
    cout<<"The bit check passed!!!"<<endl;

    daBits_phase_one.stop();
    cout << "Time required for the first phase of daBits generation and verification(for public key generation): " << daBits_phase_one.elapsed() << " seconds" << endl;

    //*RC*// cout << "used triple: " << online_op.UT.UsedTriples << endl;
    //*RC*// cout << "used square: " << online_op.UT.UsedSquares << endl;
    //*RC*// cout << "used bit: " << online_op.UT.UsedBit << endl;
    //*RC*// cout << "used input mask: " << online_op.UT.UsedInputMask << endl;
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

    // cout<<"before randbit one"<<endl;
    // P.comm_stats.print();
    for (int j = 0; j < 2 * PSIZE + SEC; j++)
    {
        for (int i = 0; i < sbit.size(); i++)
        {
            clear tmp;
            tmp = local_bits[j];
            //sbit[i].set_player(P.whoami());
            //if (i == P.whoami()) {
            online_op.get_inputs(i, sbit[i], tmp);
            
            // cout<<"after get input"<<endl;
            // P.comm_stats.print();                          
            //}
        }
        online_op.KXOR(shared_bits[j], sbit, sbit.size()); 

        // cout<<"after xor"<<endl;
        // P.comm_stats.print();  
        //shared_bits_tmp[j]=shared_bits[j];                        
    }

    // cout<<"after randbit one"<<endl;
    // P.comm_stats.print();

    //consistency check
    //r
    vector<T> r_bits(SEC+2);
    vector<T> tmp1(1);

    // cout<<"before get one"<<endl;
    // P.comm_stats.print();
    for (int i = 0; i < SEC+2; i++)
    {
      preprocessing.get_one(DATA_BIT, tmp1[0]);
      // getTuples(tmp, BIT);
      r_bits[i] = tmp1[0];
    }

    // cout<<"after get one"<<endl;
    // P.comm_stats.print();


    // vector<uint8_t> rnd1(96); // 2*381/8
    // //*RC*// P.G.get_random_bytes(rnd);
    // PRNG prng1;
    // prng1.ReSeed();
    // prng1.get_octets((octet *)rnd1.data(), sizeof(uint8_t) * rnd1.size());
    
    
    
    
    //sigma 
    
    // sigma_bits.resize(2 * PSIZE);
    // vector<T> sigma_b(2 * PSIZE);
    // vector<T> tmp_sigma(1);

    // for (int i = 0; i < sigma_b.size(); i++)
    // {
    //     preprocessing.get_one(DATA_BIT, tmp_sigma[0]);
    //     sigma_b[i]=tmp_sigma[0];
    //     //sigma_bits[i] =(rnd1[i / 8] >> (i % 8)) & 1  ;
    // }
    // vector<clear> sigma_reveal;
    // online_op.reveal({sigma_b},sigma_reveal);
    // for (int i = 0; i < sigma_bits.size(); i++)
    // {
    //     to_bigint(sigma_bits[i], sigma_reveal[i]);
    // }

    //sigma 
    vector<bigint> tmp_sigma(2 * PSIZE);
    vector<uint8_t> rnd1(96); // 2*381/8
    //*RC*// P.G.get_random_bytes(rnd);
    PRNG prng1;
    prng1.ReSeed();
    prng1.get_octets((octet *)rnd1.data(), sizeof(uint8_t) * rnd1.size());
    for (int i = 0; i < tmp_sigma.size(); i++)
    {

        tmp_sigma[i] =(rnd1[i / 8] >> (i % 8)) & 1  ;
    }
    stringstream ss;
    for(int i=0; i < tmp_sigma.size(); i++){
        ss << tmp_sigma[i];
    }
    string s(ss.str());
    // for(int i=0;i<s.length();i++){
    //     cout<<s[i];
    // }

    octetStream os_comm;
    octetStream os_open;
    octetStream os_data(s.size(), (const octet *)s.data());
    Commit(os_comm, os_open, os_data, 1); 
    // send my commitment

    // cout<<"before send all"<<endl;
    // P.comm_stats.print();
    P.send_all(os_comm); 
    // cout<<"after send all"<<endl;
    // P.comm_stats.print();

    // receive other commitments
    vector<octetStream> os_comms(P.num_players()); // CommAux(P.num_players());
    P.receive_all(os_comms); 

    // cout<<"after receive all"<<endl;
    // P.comm_stats.print();

    // send my opening
    P.send_all(os_open);

    // cout<<"after send all open"<<endl;
    // P.comm_stats.print();

    // receive other opening
    vector<octetStream> os_opens(P.num_players()); // OpenAux(P.num_players());
    P.receive_all(os_opens);

    // cout<<"after receive all open"<<endl;
    // P.comm_stats.print();


    
    octetStream os_ss;
    vector<string> sigma_share(P.num_players());
    sigma_share[P.my_num()]= s;
    for(int i=0; i<P.num_players();i++){
        if(i!=P.my_num()){
            bool res=Open(os_ss, os_comms[i],os_opens[i],1);
            if(!res){
                throw bad_dabits_value();
            }
            sigma_share[i].assign((char *)os_ss.get_data(), os_ss.get_length());
        }
    }


    sigma_bits.resize(2 * PSIZE);


    // vector<clear> neww;
    // neww.resize(762);
    vector<vector<bigint> > sigma_share_tmp(P.num_players(),vector<bigint> (762));
    for(int i=0; i<762;i++){
        for(int j=0;j<P.num_players();j++){
            string s_tmp;
            s_tmp.resize(1);
            s_tmp[0]=sigma_share[j][i];
            mpz_class bnx1(s_tmp, 2);
            bigint bn1(bnx1);
            sigma_share_tmp[j][i]=bn1;
            sigma_bits[i]=sigma_bits[i] ^ sigma_share_tmp[j][i];
        }

        // cout<< sigma_bits[i];
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
    
    // cout<<"before reveal"<<endl;
    // P.comm_stats.print();
    online_op.reveal({x_bits},x_tmp);
    
    // cout<<"after reveal"<<endl;
    // P.comm_stats.print();
    for(int i=0;i<SEC;i++){
        to_bigint(x_bits_tmp[i], x_tmp[i]);
        x_bits_plain[i]=x_bits_tmp[i] % 2;
    }


    cout<<"Perform bit check of daBits"<<endl;
    
    for(int i=0;i<2 * PSIZE;i++)
    {
        clear v=1;
        T u;
        online_op.sub_plain_new(u,v,shared_bits[i]);
        online_op.mul_inplace(u,shared_bits[i]);
        clear c;
        online_op.reveal(u,c);
        if(c != 0)
        {
            throw bad_dabits_value();
        }
    }
    cout<<"The bit check passed!!!"<<endl;

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
void POC<T>::decompose_and_reveal_for_pk_new(
    vector<bigint> &reveal_bits, const T &sk_share, const vector<T> &shared_bits, 
    int online_num, Player &P, Config_Info &CI)
{

    PRINT_DEBUG_INFO();
    reveal_bits.resize(RSIZE);

    OnlineOp<T> online_op(P, protocol, preprocessing, processor, output);
    PRINT_DEBUG_INFO();

    vector<T> deco_bits;
    PRINT_DEBUG_INFO();
    // cout<<"before a2b"<<endl;
    // P.comm_stats.print();

    online_op.A2B_new(deco_bits, sk_share);

    // cout<<"after a2b1"<<endl;
    // P.comm_stats.print();

    // cout<<"before xor"<<endl;
    // P.comm_stats.print();
    for (int i = 0; i < RSIZE; i++)
    {
        online_op.XOR_inplace(deco_bits[i], shared_bits[i]);
    }
    // cout<<"after xor"<<endl;
    // P.comm_stats.print();

    PRINT_DEBUG_INFO();
    vector<clear> out;

    // cout<<"before reveal"<<endl;
    // P.comm_stats.print();
    online_op.reveal(deco_bits, out);

    // cout<<"after reveal"<<endl;
    // P.comm_stats.print();
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
    // cout<<"before a2b"<<endl;
    // P.comm_stats.print();

    online_op.A2B(deco_bits, keys[0]);

    // cout<<"after a2b1"<<endl;
    // P.comm_stats.print();

    PRINT_DEBUG_INFO();
    online_op.A2B(tmp, keys[1]);
    
    // cout<<"after a2b2"<<endl;
    // P.comm_stats.print();

    PRINT_DEBUG_INFO();
    deco_bits.insert(deco_bits.end(), tmp.begin(), tmp.end());

    PRINT_DEBUG_INFO();
    // cout<<"before xor"<<endl;
    // P.comm_stats.print();
    for (int i = 0; i < 2 * PSIZE; i++)
    {
        online_op.XOR_inplace(deco_bits[i], shared_bits[i]);
    }
    // cout<<"after xor"<<endl;
    // P.comm_stats.print();

    PRINT_DEBUG_INFO();
    vector<clear> out;

    // cout<<"before reveal"<<endl;
    // P.comm_stats.print();
    online_op.reveal(deco_bits, out);

    // cout<<"after reveal"<<endl;
    // P.comm_stats.print();
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
   

    // cout<<"before before"<<endl;
    // P.comm_stats.print();
    for (int j = 0; j < 2 * PSIZE + SEC; j++)
    {
        // cout<<"get input before"<<endl;
        // P.comm_stats.print();

        for (int i = 0; i < sbit.size(); i++)
        {
            clear tmp;
            // tmp.assign(local_bits[j]);
            tmp = local_bits[j];
            //sbit[i].set_player(P.whoami());
            //if (i == P.whoami()) {
            online_op.get_inputs(i, sbit[i], tmp);

           
            // cout<<"get input after"<<endl;
            // P.comm_stats.print();
            //}
        }
       
        // cout<<"xor before"<<endl;
        // P.comm_stats.print();
        online_op.KXOR(shared_bits[j], sbit, sbit.size());
        // cout<<"xor after"<<endl;
        // P.comm_stats.print();
    }

   

    // cout<<"before before"<<endl;
    // P.comm_stats.print();


    //consistency Check

    //t
    vector<T> t_bits(SEC+2);
    vector<T> tmp1(1);
    // cout<<"get one before"<<endl;
    // P.comm_stats.print();
    for (int i = 0; i < SEC+2; i++)
    {
      preprocessing.get_one(DATA_BIT, tmp1[0]);
      // getTuples(tmp, BIT);
      t_bits[i] = tmp1[0];
    //   cout<<"get one after"<<endl;
    //   P.comm_stats.print();
    }

    // cout<<"get one after"<<endl;
    // P.comm_stats.print();

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
    // cout<<"reveal one before"<<endl;
    // P.comm_stats.print();
    online_op.reveal({y_bits},y_tmp);
    // cout<<"reveal one after"<<endl;
    // P.comm_stats.print();

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

    cout<<"Perform bit check of daBits"<<endl;
    
    for(int i=0;i<2 * PSIZE;i++)
    {
        clear v=1;
        T u;
        online_op.sub_plain_new(u,v,shared_bits[i]);
        online_op.mul_inplace(u,shared_bits[i]);
        clear c;
        online_op.reveal(u,c);
        if(c != 0)
        {
            throw bad_dabits_value();
        }
    }
    cout<<"The bit check passed!!!"<<endl;

    daBits_phase_two.stop();
    cout << "Time required for the second phase of daBits generation and verification: " << daBits_phase_two.elapsed() << " seconds" << endl;

    //*RC*// cout << "used triple: " << online_op.UT.UsedTriples << endl;
    //*RC*// cout << "used square: " << online_op.UT.UsedSquares << endl;
    //*RC*// cout << "used bit: " << online_op.UT.UsedBit << endl;
    //*RC*// cout << "used input mask: " << online_op.UT.UsedInputMask << endl;
}


template <class T>
void POC<T>::shared_rand_bits_for_pk_phase_two_new(
    vector<T> &shared_bits, const vector<bigint> &local_bits,const vector<bigint> &sigma_bits,
    const vector<bigint> &x_bits, int online_num, Player &P,Config_Info &CI)
{
    Timer daBits_phase_two;
    daBits_phase_two.start();
    shared_bits.resize(RSIZE+SEC);
    OnlineOp<T> online_op(P, protocol, preprocessing, processor, output);

    vector<T> sbit(P.num_players());

    // cout<<"before before"<<endl;
    // P.comm_stats.print();
    for (int j = 0; j < RSIZE + SEC; j++)
    {
        // cout<<"get input before"<<endl;
        // P.comm_stats.print();

        for (int i = 0; i < sbit.size(); i++)
        {
            clear tmp;
            // tmp.assign(local_bits[j]);
            tmp = local_bits[j];
            //sbit[i].set_player(P.whoami());
            //if (i == P.whoami()) {
            online_op.get_inputs(i, sbit[i], tmp);
            // cout<<"get input after"<<endl;
            // P.comm_stats.print();
            //}
        }
        // cout<<"xor before"<<endl;
        // P.comm_stats.print();
        online_op.KXOR(shared_bits[j], sbit, sbit.size());
        // cout<<"xor after"<<endl;
        // P.comm_stats.print();
    }

    // cout<<"before before"<<endl;
    // P.comm_stats.print();


    //consistency Check

    //t
    vector<T> t_bits(SEC+2);
    vector<T> tmp1(1);
    // cout<<"get one before"<<endl;
    // P.comm_stats.print();
    for (int i = 0; i < SEC+2; i++)
    {
      preprocessing.get_one(DATA_BIT, tmp1[0]);
      // getTuples(tmp, BIT);
      t_bits[i] = tmp1[0];
    //   cout<<"get one after"<<endl;
    //   P.comm_stats.print();
    }

    // cout<<"get one after"<<endl;
    // P.comm_stats.print();

    //y
    vector<T> y_bits(SEC);
    vector<T> tmp2(RSIZE);
    vector<T> tmp3(SEC);
    bigint TWO(2);
    for(int i=0;i<SEC;i++)
    {
        for(int j=0;j< RSIZE;j++){
            online_op.mul_plain(tmp2[j],shared_bits[j],sigma_bits[j]);
            online_op.add_inplace(y_bits[i],tmp2[j]);
        }
        online_op.add_inplace(y_bits[i],shared_bits[RSIZE+i]);
        
       
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
    // cout<<"reveal one before"<<endl;
    // P.comm_stats.print();
    online_op.reveal({y_bits},y_tmp);
    // cout<<"reveal one after"<<endl;
    // P.comm_stats.print();

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

    cout<<"Perform bit check of daBits"<<endl;
    
    for(int i=0;i<RSIZE;i++)
    {
        clear v=1;
        T u;
        online_op.sub_plain_new(u,v,shared_bits[i]);
        online_op.mul_inplace(u,shared_bits[i]);
        clear c;
        online_op.reveal(u,c);
        if(c != 0)
        {
            throw bad_dabits_value();
        }
    }
    cout<<"The bit check passed!!!"<<endl;
    daBits_phase_two.stop();
    cout << "Time required for the second phase of daBits generation and verification (pk): " << daBits_phase_two.elapsed() << " seconds" << endl;

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
void POC<T>::xor_for_pk_and_ek_new(
    vector<T> &srk, const vector<T> &shared_bits, const vector<bigint> &reveal_bits,
    int online_num, Player &P, Config_Info &CI)
{
    
    OnlineOp<T> online_op(P, protocol, preprocessing, processor, output);

    srk.resize(reveal_bits.size());
    for (int i = 0; i < reveal_bits.size(); i++)
    {
        online_op.XOR_plain(srk[i], shared_bits[i], reveal_bits[i]);
    }
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

    // cout<<"sqr before1"<<endl;
    // P.comm_stats.print();
    online_op.sqr(s0_cube, keys[0]);
    // cout<<"sqr after1"<<endl;
    // P.comm_stats.print();
    online_op.mul_inplace(s0_cube, keys[0]);

    // cout<<"sqr before2"<<endl;
    // P.comm_stats.print();

    online_op.sqr(s1_cube, keys[1]);
    // cout<<"sqr after2"<<endl;
    // P.comm_stats.print();
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

    // int res1=0;

    // for (int i = 0; i < 10; i++)
    // {
    //     clear count(i);
    //     T in;
    //     online_op.add_plain(in, uhf_out, count);
    //     int out=online_op.legendre_prf(key, in);
    //     cout<<"MPC bit"<<out<<endl;
    //     res1 &= out;
    // }
    vector<T> in;
    in.resize(10);
    for (int i = 0; i < 10; i++)
    {
        clear count(i);
        online_op.add_plain(in[i], uhf_out, count);
    }
    res=online_op.legendre_prf_new(key, in);

    //*RC*// cout << "used triple: " << online_op.UT.UsedTriples << endl;
    //*RC*// cout << "used square: " << online_op.UT.UsedSquares << endl;
    //*RC*// cout << "used bit: " << online_op.UT.UsedBit << endl;
    //*RC*// cout << "used input mask: " << online_op.UT.UsedInputMask << endl;

    return res;
}


template <class T>
int POC<T>::poc_compute_custody_bit_online_2primes_test(
    const vector<clear> pre_key, const clear key, const vector<clear> &msg, int online_num, Player &P, Config_Info &CI)
{
    if (msg.size() != CHUNK_NUM)
    {
        throw bad_value();
    }
    OnlineOp<T> online_op(P, protocol, preprocessing, processor, output);

    clear uhf_out;

    uhf_out=pre_key[0]* msg[1];
    uhf_out=uhf_out + msg[0]; // m[0] + m[1]*s1

    clear tmp;

    for (int i = 2; i < msg.size(); i++)
    {
        tmp=pre_key[i - 1]* msg[i];
        uhf_out=uhf_out + tmp;
    }

    uhf_out=uhf_out + pre_key.back();

    int out = 0;
    for (int i = 0; i < 10; i++)
    {
        clear count(i);
        clear in;
        in=uhf_out+count+key;


        int res = 0;
        
        bigint bn;

        to_bigint(bn, in);
        res = mpz_legendre(bn.get_mpz_t(), clear::pr().get_mpz_t());
        res = ceil(double(res + 1) / 2);
        cout<<"plain bit"<< res <<endl;
        out &= res;

    }

    //*RC*// cout << "used triple: " << online_op.UT.UsedTriples << endl;
    //*RC*// cout << "used square: " << online_op.UT.UsedSquares << endl;
    //*RC*// cout << "used bit: " << online_op.UT.UsedBit << endl;
    //*RC*// cout << "used input mask: " << online_op.UT.UsedInputMask << endl;

    return out;
}