#pragmaonce
// functions.h
#ifndef __FUNCTIONS_H__
#define__FUNCTIONS_H__
#include"math.h"
#include<omp.h>

externconstint N;
constdoubleV_na = -115.0, V_k = 12.0, V_l = -10.613,
g_na = 120.0, g_k = 36.0, g_l = 0.3;
constdouble a = 0.5, b = 0.1, E_syn = -10.0;
constdoublea_ex = 2, b_ex = 1, E_ex = 70.0;
doubleF_v(intnum, double* g, doubleV, double* r,
	doublem, doubleh, doublen, doubleI_ext);
doubleF_m(doubleV, doublem);
doubleF_h(doubleV, doubleh);
doubleF_n(doubleV, doublen);
doubleF_r(doubleV, doubler);
doubleF_gi(doubler_pre, doubler_post, doubleg);
#endif#pragma once
