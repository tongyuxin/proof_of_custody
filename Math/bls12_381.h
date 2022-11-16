#pragma once

#ifdef __cplusplus
extern "C"{
#endif

#include "mcl/bn_c384_256.h"

#ifdef __cplusplus
}
#endif

#include "gfp.h"
#include <string>
using namespace std;
using std::string;

const static string G1_P = "1 3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507 1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569";
const static string G2_Q = "1 352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160 3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758 1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905 927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582";
const static string ModP = "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787";

void getBasePointG1(mclBnG1 &basePoint);

void mclFr_to_G1(mclBnG1 &out, const mclBnFr &in);

void print_mclBnFr(const mclBnFr &a);

void print_mclBnG1(const mclBnG1 &a);
void print_mclBnG2(const mclBnG2 &a);

void mclBnFr_to_str(string &str, const mclBnFr &a);
void mclBnFr_to_str_new(string &str, const mclBnFr &a);
void str_to_mclBnFr(mclBnFr &out, const string &str);
void str_to_mclBnFr_new(mclBnFr &out, const string &str);

void mclBnG1_to_str(string &str, const mclBnG1 &a);
void mclBnG1_to_str_new(vector<string> &str, const mclBnG1 &a);
void str_to_mclBnG1(mclBnG1 &out, const string &str);

void mclBnG2_to_str(vector<string> &out, const mclBnG2 &a);

void mclBnG2_to_str_test(string &str, const mclBnG2 &a);

// void mclBnG2_to_gfp(vector<gfp> &out, const mclBnG2 &a);
