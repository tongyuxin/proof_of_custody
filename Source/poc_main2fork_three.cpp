#define NO_MIXED_CIRCUITS

#include "Math/gfp.hpp"
#include "Machines/SPDZ.hpp"

#include "POC/RunPoc.h"
#include "poc_args.h"
#include "ipc_msg.h"
#include <unistd.h>

template <class T>
void run(const Paras &paras);
int run_stage(const Paras &paras);

// //////////////////////
static inline std::string get_nonce(int partyid)
{
    return "123456";
}

template <typename U> // U = T::clear
inline void get_msg(int partyid, vector<U> &msg)
{
    for (size_t i = 0; i < msg.size(); i++)
        msg[i] = i * 2;
}
// //////////////////////

int main(int argc, const char **argv)
{
    Paras paras;
    if (!parse_args(argc, argv, paras))
    {
        cerr << "parse args failed!" << endl;
        exit(1);
    }

    paras.batchsize = 55;
    paras.run_stage = true;
    if (!paras.run_stage)
    {
        cout << "stage all pid:" << getpid() << " Begin!" << endl;
        paras.stage = 0;
        run_stage(paras);
        return 0;
    }

    pid_t fpid = fork();
    if (fpid < 0)
    {
        cerr << "error in fork!" << endl;
        exit(1);
    }

    int loops = 5;
    int ret = -1;
    int status = -1;
    if (fpid == 0)
    {
        pid_t fpid1 = fork();
        if (fpid1 < 0)
        {
            cerr << "error in fork!" << endl;
            exit(1);
        }

        int loops1 = 5;
        int ret1 = -1;
        int status1 = -1;

        if(fpid1==0){
            cout << "stage three pid:" << getpid() << " Begin!" << endl;
            paras.stage = 3;
            paras.lgp = paras.lgp3;
            ret = run_stage(paras);
            cout << "stage three pid:" << getpid() << " End!" << endl;

        }else{
            cout << "stage one pid:" << getpid() << " Begin!" << endl;
            paras.stage = 1;
            paras.lgp = paras.lgp1;
            paras.baseport = (paras.baseport > 50000) ? (paras.baseport - 200) : (paras.baseport + 200);
            ret1 = run_stage(paras);
            cout << "stage one pid:" << getpid() << " End!" << endl;
            waitpid(fpid1, &status1, 0);
            printf("status = %d\n", WEXITSTATUS(status1));
        }
    }
    else
    {
        cout << "stage two pid:" << getpid() << " Begin!" << endl;
        paras.stage = 2;
        paras.lgp = paras.lgp2;
        paras.baseport = (paras.baseport > 50000) ? (paras.baseport - 500) : (paras.baseport + 500);
        ret = run_stage(paras);
        cout << "stage two pid:" << getpid() << " End!" << endl;

        waitpid(fpid, &status, 0);
        printf("status = %d\n", WEXITSTATUS(status));
    }
    return 0;
}

int run_stage(const Paras &paras)
{
    int prime_length = paras.lgp;           // bit length of prime
    int n_limbs = (prime_length + 63) / 64; // compute number of 64-bit words needed
    cout << "prime_length:" << prime_length << " n_limbs:" << n_limbs << endl;

    // MASCOT
    Timer timer;
    timer.start();
    PRINT_DEBUG_INFO();
    switch (n_limbs)
    {
#define CASE_RUN(N)                    \
    case N:                            \
        run<Share<gfp_<0, N>>>(paras); \
        break
        CASE_RUN(1);
        CASE_RUN(2);
        CASE_RUN(3);
        CASE_RUN(4);
        CASE_RUN(5);
        CASE_RUN(6);
        CASE_RUN(7);
        CASE_RUN(8);
#undef CASE_RUN
    default:
        cerr << "not supported! n_limbs:" << n_limbs << endl;
        exit(1);
    }
    PRINT_DEBUG_INFO();
    cout << "POC run elapsed(s):" << timer.elapsed() << endl;

    return 0;
}

template <class T>
void run(const Paras &paras)
{
    Timer timer;
    timer.start();

    typedef typename T::clear clear;
    cout << "typeid(T).name():" << TYPENAME(typeid(T).name()) << endl;
    cout << "typeid(typename T::clear).name():" << TYPENAME(typeid(typename T::clear).name()) << endl;
    cout << "typeid(clear).name():" << TYPENAME(typeid(clear).name()) << endl;

    OnlineOptions &online_opts = OnlineOptions::singleton;
    online_opts.batch_size = paras.batchsize;
    do
    {
        string strp = "";
        if (paras.lgp == 128)
            strp = "170141183460469231731687303715885907969";
        else if (paras.lgp == 256)
            //strp = "57896044618658097711785492504343953926634992332820282019728792003956566065153";
            strp = "115792089237316195423570985008687907853269984665640564039457584007913129639747";
        else if (paras.lgp == 381)
            strp = "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787";
        else if (paras.lgp == 255)
            //strp = "727693949197446988678247011108984413331337466620091750416127";
            strp = "52435875175126190479447740508185965837690552500527637822603658699938581184513";
        else
            break;

        mpz_class bnx(strp, 10);
        bigint bn(bnx);
        online_opts.prime = bn;
    } while (0);

    Names N;
    Server::start_networking(N, paras.partyid, paras.parties, paras.hostname, paras.baseport);
    Player *player = nullptr;
    if (paras.use_encryption)
        player = new CryptoPlayer(N);
    else
        player = new PlainPlayer(N); // ThreadPlayer P(N);
    Player &P = *player;

    // initialize fields
    T::clear::init_default(paras.lgp);
    T::clear::next::init_default(paras.lgp, false);
    // must initialize MAC key for security of some protocols



    typename T::mac_key_type mac_key;
    T::read_or_generate_mac_key("", P, mac_key);
    

    // global OT setup
    BaseMachine machine;
    if (T::needs_ot)
        machine.ot_setups.push_back({P});
    

    // keeps tracks of preprocessing usage (triples etc)
    DataPositions usage;
    usage.set_num_players(P.num_players());

    // output protocol
    typename T::MAC_Check output(mac_key);
    // various preprocessing
    typename T::LivePrep preprocessing(0, usage);
    SubProcessor<T> processor(output, preprocessing, P);
    // input protocol
    typename T::Input input(processor, output);
    // multiplication protocol
    typename T::Protocol protocol(P);

    cout<<"initialization phase"<<endl;
    cout << "POC init elapsed(s):" << timer.elapsed_then_reset() << endl;
    player->comm_stats.print_and_reset();
    

    {
        POC<T> poc(P, protocol, preprocessing, processor, output);
        OnlineOp<T> online_op(P, protocol, preprocessing, processor, output);
        {
            // online_op.test_mul();
        }
        BLS<T> bls(P.num_players(), P.num_players() - 1);
        RunPOC<T> runpoc(poc, P, protocol, preprocessing, processor, output);
        Config_Info CI;

        int stage = paras.stage;

        vector<bigint> local_bits_pk, reveal_bits_pk,sigma_bits_pk,x_bits_pk;


        cout<<"setup"<<endl;
        player->comm_stats.print_and_reset();
        if(stage == 3){
                    //stage 0-0
            T sk_share;
            runpoc.run_poc_setup_new(sk_share, bls, CI);

            runpoc.run_poc_compute_public_key_phase_one_new(local_bits_pk, reveal_bits_pk, bls, CI, sigma_bits_pk, x_bits_pk ,sk_share);

            cout << "POC run_poc_compute_public_key_phase_one elapsed(s):" << timer.elapsed_then_reset() << endl;
            cout<<"public_key_phase_one"<<endl;
            player->comm_stats.print_and_reset();
            PRINT_DEBUG_INFO();


            stringstream ss;

            for (size_t i = 0; i < local_bits_pk.size(); i++)
            {
                ss << local_bits_pk[i];
            }

            for (size_t i = 0; i < reveal_bits_pk.size(); i++)
            {
                ss << reveal_bits_pk[i];
            }

            for (size_t i = 0; i < sigma_bits_pk.size(); i++)
            {
                ss << sigma_bits_pk[i];

            }

            for (size_t i = 0; i < x_bits_pk.size(); i++)
            {
                ss << x_bits_pk[i];
            }


            int size = local_bits_pk.size();
            string s(ss.str());
            int len = s.size();

            cout << "send len:" << len << ",size:" << size << endl;

            ipcmsg msg;
            msg.pathname = string("/tmp/ipc.pocnew.p" + to_string(P.my_num()));
            ipc_msg_init(msg);

            msg.buf.mtype = 1; // vector size
            msg.buf.mlen = sizeof(int);
            memcpy(msg.buf.mtext, (char *)&size, sizeof(int));
            ipc_msg_send(msg);

            msg.buf.mtype = 2; // content size
            msg.buf.mlen = sizeof(int);
            memcpy(msg.buf.mtext, (char *)&len, sizeof(int));
            ipc_msg_send(msg);

            msg.buf.mtype = 3; // content
            msg.buf.mlen = len;
            memcpy(msg.buf.mtext, s.data(), len);
            ipc_msg_send(msg);

            msg.buf.mtype = 4; // end
            ipc_msg_recv(msg);
            ipc_msg_uninit(msg);

        }

        vector<bigint> local_bits, reveal_bits,sigma_bits,x_bits;
        if (stage == 0 || stage == 1)
        {
            // string nonce = get_nonce(P.my_num());

            // stage 1-1
            //runpoc.run_poc_compute_ephem_key_2primes_phase_one(local_bits, reveal_bits, bls, nonce, CI);
            
            
            // runpoc.run_poc_compute_ephem_key_2primes_phase_one_new(local_bits, reveal_bits, bls, nonce, CI, sigma_bits, x_bits);
            // cout << "POC run_poc_compute_ephem_key_2primes_phase_one elapsed(s):" << timer.elapsed_then_reset() << endl;
            // PRINT_DEBUG_INFO();

            // cout<<"ephem_key_phase_one"<<endl;
            // player->comm_stats.print_and_reset();

            if (stage == 1)
            {

                int size;
                int len;
                string s;

                ipcmsg msg;
                msg.pathname = string("/tmp/ipc.pocnew.p" + to_string(P.my_num()));
                ipc_msg_init(msg);

                msg.buf.mtype = 1; // vector size
                ipc_msg_recv(msg);
                memcpy((char *)&size, msg.buf.mtext, msg.buf.mlen);

                msg.buf.mtype = 2; // content size
                ipc_msg_recv(msg);
                memcpy((char *)&len, msg.buf.mtext, msg.buf.mlen);

                s.resize(len);
                msg.buf.mtype = 3; // content
                ipc_msg_recv(msg);
                memcpy((char *)s.data(), msg.buf.mtext, msg.buf.mlen);
                cout << "recv len:" << len << ",size:" << size << endl;

                //int size = 295;
                local_bits_pk.resize(size);
                //reveal_bits.resize(size);
                reveal_bits_pk.resize(size - SEC);
                sigma_bits_pk.resize(RSIZE);
                x_bits_pk.resize(SEC);

                
                for (size_t i = 0; i < local_bits_pk.size(); i++)
                {

                    string s_tmp;
                    s_tmp.resize(1);
                    s_tmp[0]=s[i];
                    //s_tmp[0]='1';
                    mpz_class bnx1(s_tmp, 2);
                    bigint bn1(bnx1);
                    local_bits_pk[i]=bn1;
                    
                }
                for (size_t i = 0; i < reveal_bits_pk.size(); i++)
                {
                    string s_tmp;
                    s_tmp.resize(1);
                    s_tmp[0]=s[i+local_bits_pk.size()];
                    //s_tmp[0]='1';
                    mpz_class bnx1(s_tmp, 2);
                    bigint bn1(bnx1);
                    reveal_bits_pk[i]=bn1;
                    
                }

                for (size_t i = 0; i < sigma_bits_pk.size(); i++)
                {

                    string s_tmp;
                    s_tmp.resize(1);
                    s_tmp[0]=s[i+local_bits_pk.size()+reveal_bits_pk.size()];
                    //s_tmp[0]='1';
                    mpz_class bnx1(s_tmp, 2);
                    bigint bn1(bnx1);
                    sigma_bits_pk[i]=bn1;
                }


                for (size_t i = 0; i < x_bits_pk.size(); i++)
                {
                    string s_tmp;
                    s_tmp.resize(1);
                    s_tmp[0]=s[i+local_bits_pk.size()+reveal_bits_pk.size()+sigma_bits_pk.size()];
                    //s_tmp[0]='1';
                    mpz_class bnx1(s_tmp, 2);
                    bigint bn1(bnx1);
                    x_bits_pk[i]=bn1;
                }

                // // // separation ..................

                msg.buf.mtype = 4; // end
                msg.buf.mlen = 1;
                ipc_msg_send(msg);

                string nonce = get_nonce(P.my_num());
                vector<T> ek_tmp(2);
                runpoc.run_poc_compute_public_key_and_ek_phase_two_new(ek_tmp, bls, local_bits_pk,sigma_bits_pk,x_bits_pk, reveal_bits_pk, CI,nonce);
                cout<<"compute public key and ephem_key_phase_two"<< timer.elapsed_then_reset() << endl;
                player->comm_stats.print_and_reset();
                PRINT_DEBUG_INFO();

                runpoc.run_poc_compute_ephem_key_2primes_phase_one_new1(ek_tmp,local_bits, reveal_bits, bls, nonce, CI, sigma_bits, x_bits);
                cout << "POC run_poc_compute_ephem_key_2primes_phase_one elapsed(s):" << timer.elapsed_then_reset() << endl;
                
                cout<<"ephem_key_phase_one"<<endl;
                player->comm_stats.print_and_reset();
                // PRINT_DEBUG_INFO();


                // // separation ..................
                stringstream ss;
                for (size_t i = 0; i < local_bits.size(); i++)
                {
                    ss << local_bits[i];
                }
                for (size_t i = 0; i < reveal_bits.size(); i++)
                {
                    ss << reveal_bits[i];
                }
                for (size_t i = 0; i < sigma_bits.size(); i++)
                {
                    ss << sigma_bits[i];
                }
                for (size_t i = 0; i < x_bits.size(); i++)
                {
                    ss << x_bits[i];
                }

                int size1 = local_bits.size();
                string s1(ss.str());
                int len1 = s1.size();

                cout << "send len:" << len1 << ",size:" << size1 << endl;

                ipcmsg msg1;
                msg1.pathname = string("/tmp/ipc.poc.p" + to_string(P.my_num()));
                ipc_msg_init(msg1);

                msg1.buf.mtype = 1; // vector size
                msg1.buf.mlen = sizeof(int);
                memcpy(msg1.buf.mtext, (char *)&size1, sizeof(int));
                ipc_msg_send(msg1);

                msg1.buf.mtype = 2; // content size
                msg1.buf.mlen = sizeof(int);
                memcpy(msg1.buf.mtext, (char *)&len1, sizeof(int));
                ipc_msg_send(msg1);

                msg1.buf.mtype = 3; // content
                msg1.buf.mlen = len1;
                memcpy(msg1.buf.mtext, s1.data(), len1);
                ipc_msg_send(msg1);

                msg1.buf.mtype = 4; // end
                ipc_msg_recv(msg1);
                ipc_msg_uninit(msg1);
            }
         }
         if (stage == 0 || stage == 2)
         {
            if (stage == 2)
            {
                // cout << "stage2222222222222"  << endl;
                int size;
                int len;
                string s;

                ipcmsg msg;
                msg.pathname = string("/tmp/ipc.poc.p" + to_string(P.my_num()));
                ipc_msg_init(msg);

                msg.buf.mtype = 1; // vector size
                ipc_msg_recv(msg);
                memcpy((char *)&size, msg.buf.mtext, msg.buf.mlen);

                msg.buf.mtype = 2; // content size
                ipc_msg_recv(msg);
                memcpy((char *)&len, msg.buf.mtext, msg.buf.mlen);

                s.resize(len);
                msg.buf.mtype = 3; // content
                ipc_msg_recv(msg);
                memcpy((char *)s.data(), msg.buf.mtext, msg.buf.mlen);
                cout << "recv len:" << len << ",size:" << size << endl;

                local_bits.resize(size);
                //reveal_bits.resize(size);
                reveal_bits.resize(size - SEC);
                sigma_bits.resize(2*PSIZE);
                x_bits.resize(SEC);


                for (size_t i = 0; i < local_bits.size(); i++)
                {

                    string s_tmp;
                    s_tmp.resize(1);
                    s_tmp[0]=s[i];
                    mpz_class bnx1(s_tmp, 2);
                    bigint bn1(bnx1);
                    local_bits[i]=bn1;
                }

                for (size_t i = 0; i < reveal_bits.size(); i++)
                {
                    string s_tmp;
                    s_tmp.resize(1);
                    s_tmp[0]=s[i+local_bits.size()];
                    mpz_class bnx1(s_tmp, 2);
                    bigint bn1(bnx1);
                    reveal_bits[i]=bn1;

                }


                for (size_t i = 0; i < sigma_bits.size(); i++)
                {

                    string s_tmp;
                    s_tmp.resize(1);
                    s_tmp[0]=s[i+local_bits.size()+reveal_bits.size()];
                    mpz_class bnx1(s_tmp, 2);
                    bigint bn1(bnx1);
                    sigma_bits[i]=bn1;

                }



                for (size_t i = 0; i < x_bits.size(); i++)
                {
                    string s_tmp;
                    s_tmp.resize(1);
                    s_tmp[0]=s[i+local_bits.size()+reveal_bits.size()+sigma_bits.size()];
                    mpz_class bnx1(s_tmp, 2);
                    bigint bn1(bnx1);
                    x_bits[i]=bn1;

                }

                // // separation ..................

                msg.buf.mtype = 4; // end
                msg.buf.mlen = 1;
                ipc_msg_send(msg);
            }

        //     // online_opts.batch_size=12000;
        //     // cout<<"mul before"<<endl;
        //     // player->comm_stats.print_and_reset();
        //     // T a,b,c;
        //     // online_op.mul(a,b,c);
        //     // cout<<"mul after"<<endl;
        //     // player->comm_stats.print_and_reset();




            // stage 2-1
            vector<T> ek(3);
            //runpoc.run_poc_compute_ephem_key_2primes_phase_two(ek, local_bits, reveal_bits, CI);
            runpoc.run_poc_compute_ephem_key_2primes_phase_two_new(ek, local_bits,sigma_bits,x_bits, reveal_bits, CI);
            cout << "POC run_poc_compute_ephem_key_2primes_phase_two elapsed(s):" << timer.elapsed_then_reset() << endl;
            PRINT_DEBUG_INFO();

            cout<<"ephem_key_phase_two"<<endl;
            player->comm_stats.print_and_reset();

            // stage 2-2
            vector<T> pre_key;
            runpoc.run_poc_compute_custody_bit_offline_2primes(pre_key, ek, CI);
            cout << "POC run_poc_compute_custody_bit_offline_2primes elapsed(s):" << timer.elapsed_then_reset() << endl;
            PRINT_DEBUG_INFO();

            cout<<"offline prekey"<<endl;
            player->comm_stats.print_and_reset();

            // stage 2-3
            // for (int k = 0; k < 10; k++)
            // {
                vector<clear> msg(CHUNK_NUM);
                get_msg(P.my_num(), msg);

                PRINT_DEBUG_INFO();
                int bit = runpoc.run_poc_compute_custody_bit_online_2primes(pre_key, ek[0], msg, CI);
                cout << "POC run_poc_compute_custody_bit_online_2primes elapsed(s):" << timer.elapsed_then_reset() << endl;
                cout << "custody bit: " << bit << endl;
                PRINT_DEBUG_INFO();

                cout<<"custodybit phase"<<endl;
                player->comm_stats.print_and_reset();
            // }

         }
    }

    output.Check(P);

    cout << "endend" << endl;
    //      << std::flush;
    // player->comm_stats.print();
    // cout << "--------------------------" << endl
    //      << std::flush;
    // auto xx = preprocessing.triple_generator;
    // for (auto &p : xx->players)
    // {
    //     cout << std::flush;
    //     cout << "------------------p" << endl;
    //     usleep(1000);
    //     cout << "----------sent:" << p->sent << ", my_real_num:" << p->my_real_num()
    //          << ", elapsed:" << p->timer.elapsed() << endl;
    //     p->comm_stats.print();
    //     cout << "------------------x" << endl;
    // }

    T::LivePrep::teardown();
#ifndef VERBOSE
    //player->comm_stats.print();
#endif
    delete player;
}
