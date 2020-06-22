//Author: Chris McGinn via executable (20191029)

#ifndef PDGTOCHARGEMASS_H
#define PDGTOCHARGEMASS_H

//cpp
#include <iostream>
#include <map>

class pdgToChargeMass{
 public:
  pdgToChargeMass();
  int GetChargeFromPDG(int inPDG);
  double GetMassFromPDG(int inPDG);

 private:
  std::map<int, int> pdgToChargeMap;
  std::map<int, double> pdgToMassMap;
};

pdgToChargeMass::pdgToChargeMass(){
  //Setting map for PDG -9000211
  pdgToChargeMap[-9000211] = -1;
  pdgToMassMap[-9000211] = 1.06013;

  //Setting map for PDG -20433
  pdgToChargeMap[-20433] = -1;
  pdgToMassMap[-20433] = 2.46178;

  //Setting map for PDG -20423
  pdgToChargeMap[-20423] = 0;
  pdgToMassMap[-20423] = 2.55519;

  //Setting map for PDG -20413
  pdgToChargeMap[-20413] = -1;
  pdgToMassMap[-20413] = 2.74856;

  //Setting map for PDG -20323
  pdgToChargeMap[-20323] = -1;
  pdgToMassMap[-20323] = 1.37438;

  //Setting map for PDG -20313
  pdgToChargeMap[-20313] = 0;
  pdgToMassMap[-20313] = 1.55649;

  //Setting map for PDG -20213
  pdgToChargeMap[-20213] = -1;
  pdgToMassMap[-20213] = 1.14682;

  //Setting map for PDG -14122
  pdgToChargeMap[-14122] = -1;
  pdgToMassMap[-14122] = 2.59419;

  //Setting map for PDG -10433
  pdgToChargeMap[-10433] = -1;
  pdgToMassMap[-10433] = 2.53559;

  //Setting map for PDG -10431
  pdgToChargeMap[-10431] = -1;
  pdgToMassMap[-10431] = 2.36221;

  //Setting map for PDG -10423
  pdgToChargeMap[-10423] = 0;
  pdgToMassMap[-10423] = 2.36369;

  //Setting map for PDG -10421
  pdgToChargeMap[-10421] = 0;
  pdgToMassMap[-10421] = 2.37576;

  //Setting map for PDG -10413
  pdgToChargeMap[-10413] = -1;
  pdgToMassMap[-10413] = 2.42523;

  //Setting map for PDG -10411
  pdgToChargeMap[-10411] = -1;
  pdgToMassMap[-10411] = 2.23697;

  //Setting map for PDG -10323
  pdgToChargeMap[-10323] = -1;
  pdgToMassMap[-10323] = 1.272;

  //Setting map for PDG -10321
  pdgToChargeMap[-10321] = -1;
  pdgToMassMap[-10321] = 1.03664;

  //Setting map for PDG -10313
  pdgToChargeMap[-10313] = 0;
  pdgToMassMap[-10313] = 1.272;

  //Setting map for PDG -10311
  pdgToChargeMap[-10311] = 0;
  pdgToMassMap[-10311] = 1.27146;

  //Setting map for PDG -5324
  pdgToChargeMap[-5324] = 0;
  pdgToMassMap[-5324] = 5.97;

  //Setting map for PDG -5322
  pdgToChargeMap[-5322] = 0;
  pdgToMassMap[-5322] = 5.96;

  //Setting map for PDG -5314
  pdgToChargeMap[-5314] = 1;
  pdgToMassMap[-5314] = 5.97;

  //Setting map for PDG -5312
  pdgToChargeMap[-5312] = 1;
  pdgToMassMap[-5312] = 5.96;

  //Setting map for PDG -5232
  pdgToChargeMap[-5232] = 0;
  pdgToMassMap[-5232] = 5.788;

  //Setting map for PDG -5224
  pdgToChargeMap[-5224] = -1;
  pdgToMassMap[-5224] = 5.8321;

  //Setting map for PDG -5222
  pdgToChargeMap[-5222] = -1;
  pdgToMassMap[-5222] = 5.8113;

  //Setting map for PDG -5214
  pdgToChargeMap[-5214] = 0;
  pdgToMassMap[-5214] = 5.81;

  //Setting map for PDG -5212
  pdgToChargeMap[-5212] = 0;
  pdgToMassMap[-5212] = 5.8;

  //Setting map for PDG -5132
  pdgToChargeMap[-5132] = 1;
  pdgToMassMap[-5132] = 5.7911;

  //Setting map for PDG -5122
  pdgToChargeMap[-5122] = 0;
  pdgToMassMap[-5122] = 5.6194;

  //Setting map for PDG -5114
  pdgToChargeMap[-5114] = 1;
  pdgToMassMap[-5114] = 5.8351;

  //Setting map for PDG -5112
  pdgToChargeMap[-5112] = 1;
  pdgToMassMap[-5112] = 5.8155;

  //Setting map for PDG -4334
  pdgToChargeMap[-4334] = 0;
  pdgToMassMap[-4334] = 2.7659;

  //Setting map for PDG -4332
  pdgToChargeMap[-4332] = 0;
  pdgToMassMap[-4332] = 2.6952;

  //Setting map for PDG -4324
  pdgToChargeMap[-4324] = -1;
  pdgToMassMap[-4324] = 2.64551;

  //Setting map for PDG -4322
  pdgToChargeMap[-4322] = -1;
  pdgToMassMap[-4322] = 2.5756;

  //Setting map for PDG -4314
  pdgToChargeMap[-4314] = 0;
  pdgToMassMap[-4314] = 2.64776;

  //Setting map for PDG -4312
  pdgToChargeMap[-4312] = 0;
  pdgToMassMap[-4312] = 2.5779;

  //Setting map for PDG -4232
  pdgToChargeMap[-4232] = -1;
  pdgToMassMap[-4232] = 2.4678;

  //Setting map for PDG -4224
  pdgToChargeMap[-4224] = -2;
  pdgToMassMap[-4224] = 2.52147;

  //Setting map for PDG -4222
  pdgToChargeMap[-4222] = -2;
  pdgToMassMap[-4222] = 2.4498;

  //Setting map for PDG -4214
  pdgToChargeMap[-4214] = -1;
  pdgToMassMap[-4214] = 2.48905;

  //Setting map for PDG -4212
  pdgToChargeMap[-4212] = -1;
  pdgToMassMap[-4212] = 2.45328;

  //Setting map for PDG -4132
  pdgToChargeMap[-4132] = 0;
  pdgToMassMap[-4132] = 2.47088;

  //Setting map for PDG -4124
  pdgToChargeMap[-4124] = -1;
  pdgToMassMap[-4124] = 2.6281;

  //Setting map for PDG -4122
  pdgToChargeMap[-4122] = -1;
  pdgToMassMap[-4122] = 2.28646;

  //Setting map for PDG -4114
  pdgToChargeMap[-4114] = 0;
  pdgToMassMap[-4114] = 2.52007;

  //Setting map for PDG -4112
  pdgToChargeMap[-4112] = 0;
  pdgToMassMap[-4112] = 2.45501;

  //Setting map for PDG -3334
  pdgToChargeMap[-3334] = 1;
  pdgToMassMap[-3334] = 1.67245;

  //Setting map for PDG -3324
  pdgToChargeMap[-3324] = 0;
  pdgToMassMap[-3324] = 1.52279;

  //Setting map for PDG -3322
  pdgToChargeMap[-3322] = 0;
  pdgToMassMap[-3322] = 1.31486;

  //Setting map for PDG -3314
  pdgToChargeMap[-3314] = 1;
  pdgToMassMap[-3314] = 1.53672;

  //Setting map for PDG -3312
  pdgToChargeMap[-3312] = 1;
  pdgToMassMap[-3312] = 1.32171;

  //Setting map for PDG -3224
  pdgToChargeMap[-3224] = -1;
  pdgToMassMap[-3224] = 1.35152;

  //Setting map for PDG -3222
  pdgToChargeMap[-3222] = -1;
  pdgToMassMap[-3222] = 1.18937;

  //Setting map for PDG -3214
  pdgToChargeMap[-3214] = 0;
  pdgToMassMap[-3214] = 1.40971;

  //Setting map for PDG -3212
  pdgToChargeMap[-3212] = 0;
  pdgToMassMap[-3212] = 1.19264;

  //Setting map for PDG -3122
  pdgToChargeMap[-3122] = 0;
  pdgToMassMap[-3122] = 1.11568;

  //Setting map for PDG -3114
  pdgToChargeMap[-3114] = 1;
  pdgToMassMap[-3114] = 1.36464;

  //Setting map for PDG -3112
  pdgToChargeMap[-3112] = 1;
  pdgToMassMap[-3112] = 1.19745;

  //Setting map for PDG -2224
  pdgToChargeMap[-2224] = -2;
  pdgToMassMap[-2224] = 1.22158;

  //Setting map for PDG -2214
  pdgToChargeMap[-2214] = -1;
  pdgToMassMap[-2214] = 1.20484;

  //Setting map for PDG -2212
  pdgToChargeMap[-2212] = -1;
  pdgToMassMap[-2212] = 0.93827;

  //Setting map for PDG -2114
  pdgToChargeMap[-2114] = 0;
  pdgToMassMap[-2114] = 1.20149;

  //Setting map for PDG -2112
  pdgToChargeMap[-2112] = 0;
  pdgToMassMap[-2112] = 0.93957;

  //Setting map for PDG -1114
  pdgToChargeMap[-1114] = 1;
  pdgToMassMap[-1114] = 1.51556;

  //Setting map for PDG -543
  pdgToChargeMap[-543] = -1;
  pdgToMassMap[-543] = 6.34;

  //Setting map for PDG -541
  pdgToChargeMap[-541] = -1;
  pdgToMassMap[-541] = 6.277;

  //Setting map for PDG -533
  pdgToChargeMap[-533] = 0;
  pdgToMassMap[-533] = 5.4154;

  //Setting map for PDG -531
  pdgToChargeMap[-531] = 0;
  pdgToMassMap[-531] = 5.36677;

  //Setting map for PDG -523
  pdgToChargeMap[-523] = -1;
  pdgToMassMap[-523] = 5.3252;

  //Setting map for PDG -521
  pdgToChargeMap[-521] = -1;
  pdgToMassMap[-521] = 5.27925;

  //Setting map for PDG -513
  pdgToChargeMap[-513] = 0;
  pdgToMassMap[-513] = 5.3252;

  //Setting map for PDG -511
  pdgToChargeMap[-511] = 0;
  pdgToMassMap[-511] = 5.27958;

  //Setting map for PDG -435
  pdgToChargeMap[-435] = -1;
  pdgToMassMap[-435] = 2.55527;

  //Setting map for PDG -433
  pdgToChargeMap[-433] = -1;
  pdgToMassMap[-433] = 2.1123;

  //Setting map for PDG -431
  pdgToChargeMap[-431] = -1;
  pdgToMassMap[-431] = 1.96849;

  //Setting map for PDG -425
  pdgToChargeMap[-425] = 0;
  pdgToMassMap[-425] = 2.48045;

  //Setting map for PDG -423
  pdgToChargeMap[-423] = 0;
  pdgToMassMap[-423] = 2.00698;

  //Setting map for PDG -421
  pdgToChargeMap[-421] = 0;
  pdgToMassMap[-421] = 1.86486;

  //Setting map for PDG -415
  pdgToChargeMap[-415] = -1;
  pdgToMassMap[-415] = 2.46465;

  //Setting map for PDG -413
  pdgToChargeMap[-413] = -1;
  pdgToMassMap[-413] = 2.01028;

  //Setting map for PDG -411
  pdgToChargeMap[-411] = -1;
  pdgToMassMap[-411] = 1.86962;

  //Setting map for PDG -325
  pdgToChargeMap[-325] = -1;
  pdgToMassMap[-325] = 1.52301;

  //Setting map for PDG -323
  pdgToChargeMap[-323] = -1;
  pdgToMassMap[-323] = 0.879048;

  //Setting map for PDG -321
  pdgToChargeMap[-321] = -1;
  pdgToMassMap[-321] = 0.49368;

  //Setting map for PDG -315
  pdgToChargeMap[-315] = 0;
  pdgToMassMap[-315] = 1.40952;

  //Setting map for PDG -313
  pdgToChargeMap[-313] = 0;
  pdgToMassMap[-313] = 0.902551;

  //Setting map for PDG -311
  pdgToChargeMap[-311] = 0;
  pdgToMassMap[-311] = 0.49761;

  //Setting map for PDG -215
  pdgToChargeMap[-215] = -1;
  pdgToMassMap[-215] = 1.37495;

  //Setting map for PDG -213
  pdgToChargeMap[-213] = -1;
  pdgToMassMap[-213] = 0.744705;

  //Setting map for PDG -211
  pdgToChargeMap[-211] = -1;
  pdgToMassMap[-211] = 0.13957;

  //Setting map for PDG -16
  pdgToChargeMap[-16] = 0;
  pdgToMassMap[-16] = 0;

  //Setting map for PDG -15
  pdgToChargeMap[-15] = 1;
  pdgToMassMap[-15] = 1.77682;

  //Setting map for PDG -14
  pdgToChargeMap[-14] = 0;
  pdgToMassMap[-14] = 0;

  //Setting map for PDG -13
  pdgToChargeMap[-13] = 1;
  pdgToMassMap[-13] = 0.10566;

  //Setting map for PDG -12
  pdgToChargeMap[-12] = 0;
  pdgToMassMap[-12] = 0;

  //Setting map for PDG -11
  pdgToChargeMap[-11] = 1;
  pdgToMassMap[-11] = 0.000511;

  //Setting map for PDG -5
  pdgToChargeMap[-5] = 0;
  pdgToMassMap[-5] = 4.8;

  //Setting map for PDG -4
  pdgToChargeMap[-4] = 0;
  pdgToMassMap[-4] = 1.5;

  //Setting map for PDG -3
  pdgToChargeMap[-3] = 0;
  pdgToMassMap[-3] = 0.5;

  //Setting map for PDG -2
  pdgToChargeMap[-2] = 0;
  pdgToMassMap[-2] = 0;

  //Setting map for PDG -1
  pdgToChargeMap[-1] = 0;
  pdgToMassMap[-1] = 0.33;

  //Setting map for PDG 1
  pdgToChargeMap[1] = 0;
  pdgToMassMap[1] = 0.33;

  //Setting map for PDG 2
  pdgToChargeMap[2] = 0;
  pdgToMassMap[2] = 0;

  //Setting map for PDG 3
  pdgToChargeMap[3] = 0;
  pdgToMassMap[3] = 0.5;

  //Setting map for PDG 4
  pdgToChargeMap[4] = 0;
  pdgToMassMap[4] = 1.5;

  //Setting map for PDG 5
  pdgToChargeMap[5] = 0;
  pdgToMassMap[5] = 4.8;

  //Setting map for PDG 11
  pdgToChargeMap[11] = -1;
  pdgToMassMap[11] = 0.000511;

  //Setting map for PDG 12
  pdgToChargeMap[12] = 0;
  pdgToMassMap[12] = 0;

  //Setting map for PDG 13
  pdgToChargeMap[13] = -1;
  pdgToMassMap[13] = 0.10566;

  //Setting map for PDG 14
  pdgToChargeMap[14] = 0;
  pdgToMassMap[14] = 0;

  //Setting map for PDG 15
  pdgToChargeMap[15] = -1;
  pdgToMassMap[15] = 1.77682;

  //Setting map for PDG 16
  pdgToChargeMap[16] = 0;
  pdgToMassMap[16] = 0;

  //Setting map for PDG 21
  pdgToChargeMap[21] = 0;
  pdgToMassMap[21] = 0;

  //Setting map for PDG 22
  pdgToChargeMap[22] = 0;
  pdgToMassMap[22] = 0;

  //Setting map for PDG 90
  pdgToChargeMap[90] = 0;
  pdgToMassMap[90] = 5020;

  //Setting map for PDG 111
  pdgToChargeMap[111] = 0;
  pdgToMassMap[111] = 0.13498;

  //Setting map for PDG 113
  pdgToChargeMap[113] = 0;
  pdgToMassMap[113] = 0.870558;

  //Setting map for PDG 130
  pdgToChargeMap[130] = 0;
  pdgToMassMap[130] = 0.49761;

  //Setting map for PDG 211
  pdgToChargeMap[211] = 1;
  pdgToMassMap[211] = 0.13957;

  //Setting map for PDG 213
  pdgToChargeMap[213] = 1;
  pdgToMassMap[213] = 0.725485;

  //Setting map for PDG 215
  pdgToChargeMap[215] = 1;
  pdgToMassMap[215] = 1.33319;

  //Setting map for PDG 221
  pdgToChargeMap[221] = 0;
  pdgToMassMap[221] = 0.54785;

  //Setting map for PDG 223
  pdgToChargeMap[223] = 0;
  pdgToMassMap[223] = 0.605806;

  //Setting map for PDG 225
  pdgToChargeMap[225] = 0;
  pdgToMassMap[225] = 1.39404;

  //Setting map for PDG 310
  pdgToChargeMap[310] = 0;
  pdgToMassMap[310] = 0.49761;

  //Setting map for PDG 311
  pdgToChargeMap[311] = 0;
  pdgToMassMap[311] = 0.49761;

  //Setting map for PDG 313
  pdgToChargeMap[313] = 0;
  pdgToMassMap[313] = 0.881808;

  //Setting map for PDG 315
  pdgToChargeMap[315] = 0;
  pdgToMassMap[315] = 1.34449;

  //Setting map for PDG 321
  pdgToChargeMap[321] = 1;
  pdgToMassMap[321] = 0.49368;

  //Setting map for PDG 323
  pdgToChargeMap[323] = 1;
  pdgToMassMap[323] = 0.916706;

  //Setting map for PDG 325
  pdgToChargeMap[325] = 1;
  pdgToMassMap[325] = 1.38572;

  //Setting map for PDG 331
  pdgToChargeMap[331] = 0;
  pdgToMassMap[331] = 0.957982;

  //Setting map for PDG 333
  pdgToChargeMap[333] = 0;
  pdgToMassMap[333] = 1.01962;

  //Setting map for PDG 411
  pdgToChargeMap[411] = 1;
  pdgToMassMap[411] = 1.86962;

  //Setting map for PDG 413
  pdgToChargeMap[413] = 1;
  pdgToMassMap[413] = 2.01028;

  //Setting map for PDG 415
  pdgToChargeMap[415] = 1;
  pdgToMassMap[415] = 2.45684;

  //Setting map for PDG 421
  pdgToChargeMap[421] = 0;
  pdgToMassMap[421] = 1.86486;

  //Setting map for PDG 423
  pdgToChargeMap[423] = 0;
  pdgToMassMap[423] = 2.00698;

  //Setting map for PDG 425
  pdgToChargeMap[425] = 0;
  pdgToMassMap[425] = 2.53751;

  //Setting map for PDG 431
  pdgToChargeMap[431] = 1;
  pdgToMassMap[431] = 1.96849;

  //Setting map for PDG 433
  pdgToChargeMap[433] = 1;
  pdgToMassMap[433] = 2.1123;

  //Setting map for PDG 435
  pdgToChargeMap[435] = 1;
  pdgToMassMap[435] = 2.60207;

  //Setting map for PDG 441
  pdgToChargeMap[441] = 0;
  pdgToMassMap[441] = 2.97758;

  //Setting map for PDG 443
  pdgToChargeMap[443] = 0;
  pdgToMassMap[443] = 3.09656;

  //Setting map for PDG 445
  pdgToChargeMap[445] = 0;
  pdgToMassMap[445] = 3.55753;

  //Setting map for PDG 511
  pdgToChargeMap[511] = 0;
  pdgToMassMap[511] = 5.27958;

  //Setting map for PDG 513
  pdgToChargeMap[513] = 0;
  pdgToMassMap[513] = 5.3252;

  //Setting map for PDG 521
  pdgToChargeMap[521] = 1;
  pdgToMassMap[521] = 5.27925;

  //Setting map for PDG 523
  pdgToChargeMap[523] = 1;
  pdgToMassMap[523] = 5.3252;

  //Setting map for PDG 531
  pdgToChargeMap[531] = 0;
  pdgToMassMap[531] = 5.36677;

  //Setting map for PDG 533
  pdgToChargeMap[533] = 0;
  pdgToMassMap[533] = 5.4154;

  //Setting map for PDG 541
  pdgToChargeMap[541] = 1;
  pdgToMassMap[541] = 6.277;

  //Setting map for PDG 543
  pdgToChargeMap[543] = 1;
  pdgToMassMap[543] = 6.34;

  //Setting map for PDG 553
  pdgToChargeMap[553] = 0;
  pdgToMassMap[553] = 9.46029;

  //Setting map for PDG 555
  pdgToChargeMap[555] = 0;
  pdgToMassMap[555] = 9.9122;

  //Setting map for PDG 1103
  pdgToChargeMap[1103] = 0;
  pdgToMassMap[1103] = 1.21875;

  //Setting map for PDG 1114
  pdgToChargeMap[1114] = -1;
  pdgToMassMap[1114] = 1.12686;

  //Setting map for PDG 2101
  pdgToChargeMap[2101] = 0;
  pdgToMassMap[2101] = 0.57933;

  //Setting map for PDG 2103
  pdgToChargeMap[2103] = 0;
  pdgToMassMap[2103] = 0.77133;

  //Setting map for PDG 2112
  pdgToChargeMap[2112] = 0;
  pdgToMassMap[2112] = 0.93957;

  //Setting map for PDG 2114
  pdgToChargeMap[2114] = 0;
  pdgToMassMap[2114] = 1.27464;

  //Setting map for PDG 2203
  pdgToChargeMap[2203] = 1;
  pdgToMassMap[2203] = 0.77133;

  //Setting map for PDG 2212
  pdgToChargeMap[2212] = 1;
  pdgToMassMap[2212] = 0.93827;

  //Setting map for PDG 2214
  pdgToChargeMap[2214] = 1;
  pdgToMassMap[2214] = 1.40538;

  //Setting map for PDG 2224
  pdgToChargeMap[2224] = 2;
  pdgToMassMap[2224] = 1.22134;

  //Setting map for PDG 3101
  pdgToChargeMap[3101] = 0;
  pdgToMassMap[3101] = 1.2698;

  //Setting map for PDG 3103
  pdgToChargeMap[3103] = 0;
  pdgToMassMap[3103] = 1.72295;

  //Setting map for PDG 3112
  pdgToChargeMap[3112] = -1;
  pdgToMassMap[3112] = 1.19745;

  //Setting map for PDG 3114
  pdgToChargeMap[3114] = -1;
  pdgToMassMap[3114] = 1.36441;

  //Setting map for PDG 3122
  pdgToChargeMap[3122] = 0;
  pdgToMassMap[3122] = 1.11568;

  //Setting map for PDG 3201
  pdgToChargeMap[3201] = 0;
  pdgToMassMap[3201] = 1.39874;

  //Setting map for PDG 3203
  pdgToChargeMap[3203] = 0;
  pdgToMassMap[3203] = 1.04325;

  //Setting map for PDG 3212
  pdgToChargeMap[3212] = 0;
  pdgToMassMap[3212] = 1.19264;

  //Setting map for PDG 3214
  pdgToChargeMap[3214] = 0;
  pdgToMassMap[3214] = 1.42668;

  //Setting map for PDG 3222
  pdgToChargeMap[3222] = 1;
  pdgToMassMap[3222] = 1.18937;

  //Setting map for PDG 3224
  pdgToChargeMap[3224] = 1;
  pdgToMassMap[3224] = 1.35824;

  //Setting map for PDG 3303
  pdgToChargeMap[3303] = 0;
  pdgToMassMap[3303] = 1.16458;

  //Setting map for PDG 3312
  pdgToChargeMap[3312] = -1;
  pdgToMassMap[3312] = 1.32171;

  //Setting map for PDG 3314
  pdgToChargeMap[3314] = -1;
  pdgToMassMap[3314] = 1.5282;

  //Setting map for PDG 3322
  pdgToChargeMap[3322] = 0;
  pdgToMassMap[3322] = 1.31486;

  //Setting map for PDG 3324
  pdgToChargeMap[3324] = 0;
  pdgToMassMap[3324] = 1.52419;

  //Setting map for PDG 3334
  pdgToChargeMap[3334] = -1;
  pdgToMassMap[3334] = 1.67245;

  //Setting map for PDG 4101
  pdgToChargeMap[4101] = 0;
  pdgToMassMap[4101] = 1.91379;

  //Setting map for PDG 4103
  pdgToChargeMap[4103] = 0;
  pdgToMassMap[4103] = 2.73208;

  //Setting map for PDG 4112
  pdgToChargeMap[4112] = 0;
  pdgToMassMap[4112] = 2.44439;

  //Setting map for PDG 4114
  pdgToChargeMap[4114] = 0;
  pdgToMassMap[4114] = 2.55894;

  //Setting map for PDG 4122
  pdgToChargeMap[4122] = 1;
  pdgToMassMap[4122] = 2.28646;

  //Setting map for PDG 4124
  pdgToChargeMap[4124] = 1;
  pdgToMassMap[4124] = 2.6281;

  //Setting map for PDG 4132
  pdgToChargeMap[4132] = 0;
  pdgToMassMap[4132] = 2.47088;

  //Setting map for PDG 4201
  pdgToChargeMap[4201] = 1;
  pdgToMassMap[4201] = 2.63104;

  //Setting map for PDG 4203
  pdgToChargeMap[4203] = 1;
  pdgToMassMap[4203] = 2.74548;

  //Setting map for PDG 4212
  pdgToChargeMap[4212] = 1;
  pdgToMassMap[4212] = 2.45393;

  //Setting map for PDG 4214
  pdgToChargeMap[4214] = 1;
  pdgToMassMap[4214] = 2.51449;

  //Setting map for PDG 4222
  pdgToChargeMap[4222] = 2;
  pdgToMassMap[4222] = 2.45434;

  //Setting map for PDG 4224
  pdgToChargeMap[4224] = 2;
  pdgToMassMap[4224] = 2.52191;

  //Setting map for PDG 4232
  pdgToChargeMap[4232] = 1;
  pdgToMassMap[4232] = 2.4678;

  //Setting map for PDG 4303
  pdgToChargeMap[4303] = 0;
  pdgToMassMap[4303] = 2.74153;

  //Setting map for PDG 4312
  pdgToChargeMap[4312] = 0;
  pdgToMassMap[4312] = 2.5779;

  //Setting map for PDG 4314
  pdgToChargeMap[4314] = 0;
  pdgToMassMap[4314] = 2.65513;

  //Setting map for PDG 4322
  pdgToChargeMap[4322] = 1;
  pdgToMassMap[4322] = 2.5756;

  //Setting map for PDG 4324
  pdgToChargeMap[4324] = 1;
  pdgToMassMap[4324] = 2.64389;

  //Setting map for PDG 4332
  pdgToChargeMap[4332] = 0;
  pdgToMassMap[4332] = 2.6952;

  //Setting map for PDG 4334
  pdgToChargeMap[4334] = 0;
  pdgToMassMap[4334] = 2.7659;

  //Setting map for PDG 4403
  pdgToChargeMap[4403] = 1;
  pdgToMassMap[4403] = 3.77632;

  //Setting map for PDG 5101
  pdgToChargeMap[5101] = 0;
  pdgToMassMap[5101] = 5.69972;

  //Setting map for PDG 5103
  pdgToChargeMap[5103] = 0;
  pdgToMassMap[5103] = 5.91342;

  //Setting map for PDG 5112
  pdgToChargeMap[5112] = -1;
  pdgToMassMap[5112] = 5.8155;

  //Setting map for PDG 5114
  pdgToChargeMap[5114] = -1;
  pdgToMassMap[5114] = 5.8351;

  //Setting map for PDG 5122
  pdgToChargeMap[5122] = 0;
  pdgToMassMap[5122] = 5.6194;

  //Setting map for PDG 5132
  pdgToChargeMap[5132] = -1;
  pdgToMassMap[5132] = 5.7911;

  //Setting map for PDG 5203
  pdgToChargeMap[5203] = 0;
  pdgToMassMap[5203] = 5.56413;

  //Setting map for PDG 5212
  pdgToChargeMap[5212] = 0;
  pdgToMassMap[5212] = 5.8;

  //Setting map for PDG 5214
  pdgToChargeMap[5214] = 0;
  pdgToMassMap[5214] = 5.81;

  //Setting map for PDG 5222
  pdgToChargeMap[5222] = 1;
  pdgToMassMap[5222] = 5.8113;

  //Setting map for PDG 5224
  pdgToChargeMap[5224] = 1;
  pdgToMassMap[5224] = 5.8321;

  //Setting map for PDG 5232
  pdgToChargeMap[5232] = 0;
  pdgToMassMap[5232] = 5.788;

  //Setting map for PDG 5314
  pdgToChargeMap[5314] = -1;
  pdgToMassMap[5314] = 5.97;

  //Setting map for PDG 5322
  pdgToChargeMap[5322] = 0;
  pdgToMassMap[5322] = 5.96;

  //Setting map for PDG 5324
  pdgToChargeMap[5324] = 0;
  pdgToMassMap[5324] = 5.97;

  //Setting map for PDG 10113
  pdgToChargeMap[10113] = 0;
  pdgToMassMap[10113] = 1.56885;

  //Setting map for PDG 10213
  pdgToChargeMap[10213] = 1;
  pdgToMassMap[10213] = 1.15707;

  //Setting map for PDG 10221
  pdgToChargeMap[10221] = 0;
  pdgToMassMap[10221] = 1.22385;

  //Setting map for PDG 10311
  pdgToChargeMap[10311] = 0;
  pdgToMassMap[10311] = 1.53444;

  //Setting map for PDG 10313
  pdgToChargeMap[10313] = 0;
  pdgToMassMap[10313] = 1.272;

  //Setting map for PDG 10321
  pdgToChargeMap[10321] = 1;
  pdgToMassMap[10321] = 1.08179;

  //Setting map for PDG 10323
  pdgToChargeMap[10323] = 1;
  pdgToMassMap[10323] = 1.272;

  //Setting map for PDG 10411
  pdgToChargeMap[10411] = 1;
  pdgToMassMap[10411] = 2.23757;

  //Setting map for PDG 10413
  pdgToChargeMap[10413] = 1;
  pdgToMassMap[10413] = 2.42538;

  //Setting map for PDG 10421
  pdgToChargeMap[10421] = 0;
  pdgToMassMap[10421] = 2.32986;

  //Setting map for PDG 10423
  pdgToChargeMap[10423] = 0;
  pdgToMassMap[10423] = 2.4297;

  //Setting map for PDG 10431
  pdgToChargeMap[10431] = 1;
  pdgToMassMap[10431] = 2.31104;

  //Setting map for PDG 10433
  pdgToChargeMap[10433] = 1;
  pdgToMassMap[10433] = 2.53486;

  //Setting map for PDG 10441
  pdgToChargeMap[10441] = 0;
  pdgToMassMap[10441] = 3.41154;

  //Setting map for PDG 10443
  pdgToChargeMap[10443] = 0;
  pdgToMassMap[10443] = 3.50165;

  //Setting map for PDG 10551
  pdgToChargeMap[10551] = 0;
  pdgToMassMap[10551] = 9.8594;

  //Setting map for PDG 14122
  pdgToChargeMap[14122] = 1;
  pdgToMassMap[14122] = 2.59044;

  //Setting map for PDG 20113
  pdgToChargeMap[20113] = 0;
  pdgToMassMap[20113] = 1.70564;

  //Setting map for PDG 20213
  pdgToChargeMap[20213] = 1;
  pdgToMassMap[20213] = 1.36702;

  //Setting map for PDG 20313
  pdgToChargeMap[20313] = 0;
  pdgToMassMap[20313] = 1.60978;

  //Setting map for PDG 20413
  pdgToChargeMap[20413] = 1;
  pdgToMassMap[20413] = 2.33784;

  //Setting map for PDG 20423
  pdgToChargeMap[20423] = 0;
  pdgToMassMap[20423] = 2.53267;

  //Setting map for PDG 20433
  pdgToChargeMap[20433] = 1;
  pdgToMassMap[20433] = 2.44013;

  //Setting map for PDG 20443
  pdgToChargeMap[20443] = 0;
  pdgToMassMap[20443] = 3.51084;

  //Setting map for PDG 20553
  pdgToChargeMap[20553] = 0;
  pdgToMassMap[20553] = 9.8928;

  //Setting map for PDG 30443
  pdgToChargeMap[30443] = 0;
  pdgToMassMap[30443] = 3.77587;

  //Setting map for PDG 100441
  pdgToChargeMap[100441] = 0;
  pdgToMassMap[100441] = 3.63976;

  //Setting map for PDG 100443
  pdgToChargeMap[100443] = 0;
  pdgToMassMap[100443] = 3.68626;

  //Setting map for PDG 9010221
  pdgToChargeMap[9010221] = 0;
  pdgToMassMap[9010221] = 1;

  //Setting map for PDG 9940003
  pdgToChargeMap[9940003] = 0;
  pdgToMassMap[9940003] = 3.29692;

  //Setting map for PDG 9940005
  pdgToChargeMap[9940005] = 0;
  pdgToMassMap[9940005] = 3.7562;

  //Setting map for PDG 9940011
  pdgToChargeMap[9940011] = 0;
  pdgToMassMap[9940011] = 3.61475;

  //Setting map for PDG 9940023
  pdgToChargeMap[9940023] = 0;
  pdgToMassMap[9940023] = 3.71066;

  //Setting map for PDG 9940103
  pdgToChargeMap[9940103] = 0;
  pdgToMassMap[9940103] = 3.88611;

  //Setting map for PDG 9941003
  pdgToChargeMap[9941003] = 0;
  pdgToMassMap[9941003] = 3.29692;

  //Setting map for PDG 9941103
  pdgToChargeMap[9941103] = 0;
  pdgToMassMap[9941103] = 3.88611;

  //Setting map for PDG 9942003
  pdgToChargeMap[9942003] = 0;
  pdgToMassMap[9942003] = 3.29692;

  //Setting map for PDG 9942033
  pdgToChargeMap[9942033] = 0;
  pdgToMassMap[9942033] = 3.97315;

  //Setting map for PDG 9942103
  pdgToChargeMap[9942103] = 0;
  pdgToMassMap[9942103] = 3.88611;

}
int pdgToChargeMass::GetChargeFromPDG(int inPDG){
 if(pdgToChargeMap.count(inPDG) == 0) std::cout << "PDGTOCHARGEMASS ERROR: given inPDG val '" << inPDG << "' is invalid. returning random from map" << std::endl;
 return pdgToChargeMap[inPDG];
}

double pdgToChargeMass::GetMassFromPDG(int inPDG){
 if(pdgToMassMap.count(inPDG) == 0) std::cout << "PDGTOCHARGEMASS ERROR: given inPDG val '" << inPDG << "' is invalid. returning random from map" << std::endl;
 return pdgToMassMap[inPDG];
}
#endif
