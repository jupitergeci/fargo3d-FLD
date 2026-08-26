/*
 * Semenov et al. (2003)
 * Rosseland mean opacity: nrm / h / s.
 *
 * Dense table generated offline with:
 *
 *   NRHO  = 481
 *   NTEMP = 4001
 *
 * Uniform axes:
 *
 *   -18 <= log10(rho_cgs) <= -7
 *   log10(5 K) <= log10(T) <= 4
 */

#ifndef RADIATION_OPACITY_TABLE_H
#define RADIATION_OPACITY_TABLE_H

#define SEMENOV_NRHO 481
#define SEMENOV_NTEMP 4001

#define SEMENOV_LOGRHO_MIN (-18.0)
#define SEMENOV_LOGRHO_MAX (-7.0)

#define SEMENOV_LOGT_MIN 0.6989700043360189
#define SEMENOV_LOGT_MAX 4.0

#define SEMENOV_T_MIN 5.0
#define SEMENOV_T_MAX 10000.0

#define SEMENOV_INV_DLOGRHO \
((SEMENOV_NRHO-1)/(SEMENOV_LOGRHO_MAX-SEMENOV_LOGRHO_MIN))

#define SEMENOV_INV_DLOGT \
((SEMENOV_NTEMP-1)/(SEMENOV_LOGT_MAX-SEMENOV_LOGT_MIN))

#define SEMENOV_TABLE_SIZE \
(SEMENOV_NRHO*SEMENOV_NTEMP)

#ifndef SEMENOV_OPACITY_FILE
#define SEMENOV_OPACITY_FILE \
"opacity/semenov_rosseland_nrm_h_s_481x4001.dat"
#endif

#ifndef RAD_OPACITY_MU
#define RAD_OPACITY_MU 2.34
#endif

#define RAD_OPACITY_R_GAS_CGS 8.31446261815324e7
#define RAD_OPACITY_R_MU_CGS \
(RAD_OPACITY_R_GAS_CGS/RAD_OPACITY_MU)

#endif