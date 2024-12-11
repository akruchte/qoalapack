/* Symbol registration initialization: original provided by Brian Ripley.
   Anything called from R should be registered here (and declared in mgcv.h).
   (See also NAMESPACE:1)
 */ 
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "qmgcv.h"

R_CallMethodDef CallMethods[] = {
  {"qmgcv_pmmult2", (DL_FUNC) &qmgcv_pmmult2,5},
  {"qmgcv_Rpiqr", (DL_FUNC) &qmgcv_Rpiqr,5}, 
  { "qmgcv_tmm",(DL_FUNC)&qmgcv_tmm,5},
  { "qmgcv_chol_down",(DL_FUNC)&qmgcv_chol_down,5},
  { "qmgcv_chol_up",(DL_FUNC)&qmgcv_chol_up,5}, 
  { "qmgcv_Rpbsi",(DL_FUNC)&qmgcv_Rpbsi,2},
  { "qmgcv_RPPt",(DL_FUNC)&qmgcv_RPPt,3},
  { "qmgcv_Rpchol",(DL_FUNC)&qmgcv_Rpchol,4},
  { "qmgcv_Rpforwardsolve",(DL_FUNC)&qmgcv_Rpforwardsolve,3},
  { "qmgcv_Rpbacksolve",(DL_FUNC)&qmgcv_Rpbacksolve,3},
  { "qmgcv_Rpcross",(DL_FUNC)&qmgcv_Rpcross,3},
  { "qmgcv_madi",(DL_FUNC)&qmgcv_madi,4},
  { "Rkdtree",(DL_FUNC)&Rkdtree,1},
  {"Rkdnearest",(DL_FUNC)&Rkdnearest,4},
  {"Rkradius",(DL_FUNC)&Rkradius,5},
  {"sXWXd",(DL_FUNC)&sXWXd,5},
  {"CXWXd0",(DL_FUNC)&CXWXd0,14},
  {"CXWXd1",(DL_FUNC)&CXWXd1,16},
  {"sXbd",(DL_FUNC)&sXbd,3},
  {"CXbd",(DL_FUNC)&CXbd,13},
  {"sXyd",(DL_FUNC)&sXyd,3},
  {"CXWyd",(DL_FUNC)&CXWyd,17},
  {"sdiagXVXt",(DL_FUNC)&sdiagXVXt,4},
  {"CdiagXVXt",(DL_FUNC)&CdiagXVXt,14},
  {"stmm",(DL_FUNC)&stmm,1},
  {"AddBVB",(DL_FUNC)&AddBVB,3},
  {"isa1p",(DL_FUNC)&isa1p,3},
  {"mrow_sum",(DL_FUNC)&mrow_sum,3},
  {"ncv",(DL_FUNC)&ncv,17},
  {"Rncv",(DL_FUNC)&Rncv,19},
  {"ncvls",(DL_FUNC)&ncvls,18},
  {"Rncvls",(DL_FUNC)&Rncvls,19},
  {"nei_cov",(DL_FUNC)&nei_cov,4},
  {"dpdev",(DL_FUNC)&dpdev,1},
  {"spdev",(DL_FUNC)&spdev,1},
  {"mvnll",(DL_FUNC)&mvnll,13},
  {NULL, NULL, 0}
};

R_CMethodDef CEntries[] = { 
    {"band_chol",(DL_FUNC) band_chol,4},
    {"davies",(DL_FUNC) davies,10},
    {"tri_chol",(DL_FUNC) tri_chol,4},
    {"diagXVXt", (DL_FUNC) &diagXVXt,21},
    {"XWXd0", (DL_FUNC) &XWXd0,18},
    {"XWXd1", (DL_FUNC) &XWXd1,22},
    {"XWyd", (DL_FUNC) &XWyd,21},
    {"Xbd", (DL_FUNC) &Xbd,17},
    {"vcorr", (DL_FUNC) &vcorr, 5},
    {"dchol", (DL_FUNC) &dchol, 4},
    {"chol_down", (DL_FUNC) &chol_down, 5},
    {"qmgcv_omp", (DL_FUNC) &qmgcv_omp, 1},
    {"coxpred", (DL_FUNC) &coxpred, 14},
    {"coxpp", (DL_FUNC) &coxpp, 10},
    {"coxlpl", (DL_FUNC) &coxlpl, 17},
    // {"mvn_ll", (DL_FUNC) &mvn_ll,15},
    {"RMonoCon", (DL_FUNC) &RMonoCon, 7},
    {"RuniqueCombs", (DL_FUNC) &RuniqueCombs, 4},
    {"RPCLS", (DL_FUNC) &RPCLS, 13},
    {"construct_tprs", (DL_FUNC) &construct_tprs, 13},
    {"crspl", (DL_FUNC) &crspl,8},
    {"predict_tprs", (DL_FUNC) &predict_tprs, 12},
    {"MinimumSeparation", (DL_FUNC) &MinimumSeparation, 6},
    {"magic", (DL_FUNC) &magic, 19},
    {"qmgcv_mmult", (DL_FUNC) &qmgcv_mmult,8},
    {"qmgcv_pmmult", (DL_FUNC) &qmgcv_pmmult,9},
    {"gdi1",(DL_FUNC) &gdi1,49},
    {"gdi2",(DL_FUNC) &gdi2,48},
    {"R_cond",(DL_FUNC) &R_cond,5} ,
    {"pls_fit1",(DL_FUNC)&pls_fit1,14},
    {"tweedious",(DL_FUNC)&tweedious,13},
    {"tweedious2",(DL_FUNC)&tweedious2,13},
    {"psum",(DL_FUNC)&psum,4},
    {"get_detS2",(DL_FUNC)&get_detS2,12},
    {"get_stableS",(DL_FUNC)&get_stableS,14},
    {"qmgcv_tri_diag",(DL_FUNC)&qmgcv_tri_diag,3},
    {"qmgcv_td_qy",(DL_FUNC)&qmgcv_td_qy,7},
    {"qmgcv_symeig",(DL_FUNC)&qmgcv_symeig,6},
    {"qmgcv_trisymeig",(DL_FUNC)&qmgcv_trisymeig,6},
    {"read_mat",(DL_FUNC)&read_mat,4},
    {"rwMatrix",(DL_FUNC)&rwMatrix,8},
    {"in_out",(DL_FUNC)&in_out,8},
    {"Rlanczos",(DL_FUNC)&Rlanczos,8},
    {"rksos",(DL_FUNC)&rksos,3},
    {"gen_tps_poly_powers",(DL_FUNC)&gen_tps_poly_powers,4},
    {"k_nn",(DL_FUNC)&k_nn,8},
    // {"Rkdtree",(DL_FUNC)&Rkdtree,5},
    //{"Rkdnearest",(DL_FUNC)&Rkdnearest,9},
    //{"Rkradius",(DL_FUNC)&Rkradius,9},
    {"sspl_construct",(DL_FUNC)&sspl_construct,9},
    {"sspl_mapply",(DL_FUNC)&sspl_mapply,9},
    {"tri2nei",(DL_FUNC)&tri2nei,5},
    {"nei_penalty",(DL_FUNC)&nei_penalty, 10},
    {"boundary",(DL_FUNC)&boundary, 14},
    {"pde_coeffs",(DL_FUNC)&pde_coeffs, 9},
    {"gridder",(DL_FUNC)&gridder, 13},
    {"row_block_reorder",(DL_FUNC)&row_block_reorder,5},
    {"qmgcv_pqr",(DL_FUNC)&qmgcv_pqr,6},
    {"getRpqr",(DL_FUNC)&getRpqr,6},
    {"qmgcv_pqrqy",(DL_FUNC)&qmgcv_pqrqy,8},
    {"minres",(DL_FUNC)&minres,7},
    {"Zb",(DL_FUNC)&Zb,6},
    {"Ztb",(DL_FUNC)&Ztb,7},
    {NULL, NULL, 0}
};

void R_init_qmgcv(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_RegisterCCallable("qmgcv","qmgcv_pmmult2", (DL_FUNC) &qmgcv_pmmult2); // allows calling from other packages
    R_RegisterCCallable("qmgcv","pls_fit1", (DL_FUNC) &pls_fit1);
    R_RegisterCCallable("qmgcv","gdi2", (DL_FUNC) &gdi2); 
}