*----------------------------------------------------------------------
* sim5.do 
*----------------------------------------------------------------------

version 16
clear all
set more off
set seed 82005

local reps           2000
local alpha           0.05
local intercept       .10
local cv              .05
local pretrend        0
local step1           0
local Nlist            "10 15 20"
local posttrend_list  "0 2 4 8"
local Tlist           "30 60 100"
local rho1val         0.7
local scen            3

postfile sim_res          ///
    int    scenario       ///  1=mild positive  2=oscillatory  3=high persistent
    int    N_panels       ///
    int    T_periods      ///
    double posttrend      ///
    int    method         ///  1=naive (oim)  2=xtbetark2 vce(dk)  3=xtbetark2 vce(pcse)
    double reject         ///  TYPE I ERROR when posttrend==0, POWER when posttrend!=0
    double coverage       ///
    double pct_bias       ///
    double rmse           ///
    double se_ratio        ///  SE relative to method=1 (naive/oim) in the SAME replication
    double rho1_bias      ///  percentage bias in the AR(1) coefficient (rho1)
    double rho1_rmse      ///  RMSE of the AR(1) coefficient (rho1)
    double rho1_coverage  ///  95% CI coverage rate for rho1
    double rho1_se_ratio  ///  SE for rho1 relative to method=1 (naive/oim) in the SAME replication
    using "E:\ArielStuff\Simulations\sim5.dta", replace

foreach Nval of numlist `Nlist' {
foreach T of numlist `Tlist' {

    local trperiod1 = floor(`T'/2)

    local eta0 = logit(`intercept')
    local mu_pre_target = `intercept' * (1 + `pretrend'/100)
    local b_pre = (logit(`mu_pre_target') - `eta0') / `trperiod1'
    local mu_cf = invlogit(`eta0' + `b_pre'*`trperiod1')
    local mu_step_target = `mu_cf' + `intercept'*`step1'/100
    local b_step1 = logit(`mu_step_target') - logit(`mu_cf')
    local post_horizon = `T' - 1 - `trperiod1'
    local phi = (1-`intercept')/(`intercept'*`cv'^2) - 1

foreach posttrend of numlist `posttrend_list' {

    local mu_post_target = `mu_step_target' + `intercept'*`posttrend'/100
    local b_post1 = (logit(`mu_post_target') - logit(`mu_step_target')) / `post_horizon'
    local b1_true = `b_post1' - `b_pre'

    di as txt _newline ///
        "Scen=`scen' (rho1=0.7)  N=`Nval'  T=`T'  posttrend=`posttrend'"

    local sum_b_naive=0
    local sum_b_dk=0
    local sum_b_pcse=0
    local sum_sq_naive=0
    local sum_sq_dk=0
    local sum_sq_pcse=0
    local sum_se_naive=0
    local sum_se_dk=0
    local sum_se_pcse=0
    local n_rej_naive=0
    local n_rej_dk=0
    local n_rej_pcse=0
    local n_cov_naive=0
    local n_cov_dk=0
    local n_cov_pcse=0
    local sum_rho1_b_naive=0
    local sum_rho1_b_dk=0
    local sum_rho1_b_pcse=0
    local sum_rho1_sq_naive=0
    local sum_rho1_sq_dk=0
    local sum_rho1_sq_pcse=0
    local sum_rho1_se_naive=0
    local sum_rho1_se_dk=0
    local sum_rho1_se_pcse=0
    local n_rho1_cov_naive=0
    local n_rho1_cov_dk=0
    local n_rho1_cov_pcse=0
    local n_ok=0
    local n_nonstat=0

    _dots 0, title(Replications) reps(`reps')
    forvalues rep = 1/`reps' {

        quietly {
            clear
            set obs `= `Nval'*`T''
            gen int    id = ceil(_n/`T')
            gen int    t  = mod(_n-1,`T')
            sort id t
            tsset id t

            gen byte   x1  = (t >= `trperiod1')
            gen double x1t = (t - `trperiod1') * x1
            gen _x`trperiod1'   = x1
            gen _x_t`trperiod1' = x1t

            gen double xb = `eta0' + `b_pre'*t + `b_step1'*x1 + (`b_post1'-`b_pre')*x1t

            gen double y = .
            by id: replace y = rbeta(invlogit(xb)*`phi', (1-invlogit(xb))*`phi') if _n==1
            by id: replace y = rbeta( ///
                    invlogit(xb + `rho1val'*(logit(y[_n-1])-xb[_n-1]))*`phi', ///
                    (1-invlogit(xb + `rho1val'*(logit(y[_n-1])-xb[_n-1])))*`phi' ///
                ) if _n>1
        }

        capture quietly xtbetark2 y _x`trperiod1' _x_t`trperiod1', lag(1) vce(all) nolog
        if _rc {
            _dots `rep' 1
            continue
        }

        local brho1 = _b[ar:rho1]

        if abs(`brho1') >= 1 {
            local n_nonstat = `n_nonstat' + 1
            _dots `rep' 1
            continue
        }

        tempname Voim Vdk Vpcse
        matrix `Voim'  = e(V_oim)
        matrix `Vdk'   = e(V_dk)
        matrix `Vpcse' = e(V_pcse)
        local vname "y:_x_t`trperiod1'"

        local bhat   = _b[`vname']
        local snaive = sqrt(`Voim'["`vname'","`vname'"])
        local sdk    = sqrt(`Vdk'["`vname'","`vname'"])
        local spcse  = sqrt(`Vpcse'["`vname'","`vname'"])
        if missing(`bhat') | missing(`snaive') | missing(`sdk') | missing(`spcse') {
            _dots `rep' 1
            continue
        }

        local serho1_naive = sqrt(`Voim'["ar:rho1","ar:rho1"])
        local serho1_dk    = sqrt(`Vdk'["ar:rho1","ar:rho1"])
        local serho1_pcse  = sqrt(`Vpcse'["ar:rho1","ar:rho1"])
        if missing(`serho1_naive') | missing(`serho1_dk') | missing(`serho1_pcse') {
            _dots `rep' 1
            continue
        }

        local pnaive  = 2*(1-normal(abs(`bhat'/`snaive')))
        local lonaive = `bhat' - invnormal(0.975)*`snaive'
        local hinaive = `bhat' + invnormal(0.975)*`snaive'
        local pdk    = 2*(1-normal(abs(`bhat'/`sdk')))
        local lodk   = `bhat' - invnormal(0.975)*`sdk'
        local hidk   = `bhat' + invnormal(0.975)*`sdk'
        local ppcse  = 2*(1-normal(abs(`bhat'/`spcse')))
        local lopcse = `bhat' - invnormal(0.975)*`spcse'
        local hipcse = `bhat' + invnormal(0.975)*`spcse'

        local lorho1_naive = `brho1' - invnormal(0.975)*`serho1_naive'
        local hirho1_naive = `brho1' + invnormal(0.975)*`serho1_naive'
        local lorho1_dk    = `brho1' - invnormal(0.975)*`serho1_dk'
        local hirho1_dk    = `brho1' + invnormal(0.975)*`serho1_dk'
        local lorho1_pcse  = `brho1' - invnormal(0.975)*`serho1_pcse'
        local hirho1_pcse  = `brho1' + invnormal(0.975)*`serho1_pcse'

        local sum_b_naive  = `sum_b_naive'  + `bhat'
        local sum_b_dk     = `sum_b_dk'     + `bhat'
        local sum_b_pcse   = `sum_b_pcse'   + `bhat'
        local sum_sq_naive = `sum_sq_naive' + (`bhat'-`b1_true')^2
        local sum_sq_dk    = `sum_sq_dk'    + (`bhat'-`b1_true')^2
        local sum_sq_pcse  = `sum_sq_pcse'  + (`bhat'-`b1_true')^2
        local sum_se_naive = `sum_se_naive' + `snaive'
        local sum_se_dk    = `sum_se_dk'    + `sdk'
        local sum_se_pcse  = `sum_se_pcse'  + `spcse'
        if `pnaive'<`alpha'                            local n_rej_naive = `n_rej_naive'+1
        if `lonaive'<=`b1_true' & `b1_true'<=`hinaive'  local n_cov_naive = `n_cov_naive'+1
        if `pdk'<`alpha'                                local n_rej_dk    = `n_rej_dk'+1
        if `lodk'<=`b1_true'   & `b1_true'<=`hidk'      local n_cov_dk    = `n_cov_dk'+1
        if `ppcse'<`alpha'                              local n_rej_pcse  = `n_rej_pcse'+1
        if `lopcse'<=`b1_true' & `b1_true'<=`hipcse'    local n_cov_pcse  = `n_cov_pcse'+1

        local sum_rho1_b_naive  = `sum_rho1_b_naive'  + `brho1'
        local sum_rho1_sq_naive = `sum_rho1_sq_naive' + (`brho1'-`rho1val')^2
        local sum_rho1_se_naive = `sum_rho1_se_naive' + `serho1_naive'
        if `lorho1_naive'<=`rho1val' & `rho1val'<=`hirho1_naive' local n_rho1_cov_naive = `n_rho1_cov_naive'+1
        local sum_rho1_b_dk  = `sum_rho1_b_dk'  + `brho1'
        local sum_rho1_sq_dk = `sum_rho1_sq_dk' + (`brho1'-`rho1val')^2
        local sum_rho1_se_dk = `sum_rho1_se_dk' + `serho1_dk'
        if `lorho1_dk'<=`rho1val' & `rho1val'<=`hirho1_dk' local n_rho1_cov_dk = `n_rho1_cov_dk'+1
        local sum_rho1_b_pcse  = `sum_rho1_b_pcse'  + `brho1'
        local sum_rho1_sq_pcse = `sum_rho1_sq_pcse' + (`brho1'-`rho1val')^2
        local sum_rho1_se_pcse = `sum_rho1_se_pcse' + `serho1_pcse'
        if `lorho1_pcse'<=`rho1val' & `rho1val'<=`hirho1_pcse' local n_rho1_cov_pcse = `n_rho1_cov_pcse'+1

        local n_ok = `n_ok' + 1
        _dots `rep' 0
    }


    local rmse_naive = sqrt(`sum_sq_naive'/`n_ok')
    local rmse_dk    = sqrt(`sum_sq_dk'/`n_ok')
    local rmse_pcse  = sqrt(`sum_sq_pcse'/`n_ok')
    local mse_naive  = `sum_se_naive'/`n_ok'
    local mse_dk     = `sum_se_dk'/`n_ok'
    local mse_pcse   = `sum_se_pcse'/`n_ok'
    local mb_naive   = `sum_b_naive'/`n_ok'
    local mb_dk      = `sum_b_dk'/`n_ok'
    local mb_pcse    = `sum_b_pcse'/`n_ok'

    local rho1_rmse_naive = sqrt(`sum_rho1_sq_naive'/`n_ok')
    local rho1_rmse_dk    = sqrt(`sum_rho1_sq_dk'/`n_ok')
    local rho1_rmse_pcse  = sqrt(`sum_rho1_sq_pcse'/`n_ok')
    local rho1_mse_naive  = `sum_rho1_se_naive'/`n_ok'
    local rho1_mse_dk     = `sum_rho1_se_dk'/`n_ok'
    local rho1_mse_pcse   = `sum_rho1_se_pcse'/`n_ok'
    local rho1_mb_naive   = `sum_rho1_b_naive'/`n_ok'
    local rho1_mb_dk      = `sum_rho1_b_dk'/`n_ok'
    local rho1_mb_pcse    = `sum_rho1_b_pcse'/`n_ok'
    local rho1_bias_naive = 100*(`rho1_mb_naive'-`rho1val')/`rho1val'
    local rho1_bias_dk    = 100*(`rho1_mb_dk'-`rho1val')/`rho1val'
    local rho1_bias_pcse  = 100*(`rho1_mb_pcse'-`rho1val')/`rho1val'

    if `b1_true'!=0 {
        local bias_naive = 100*(`mb_naive'-`b1_true')/`b1_true'
        local bias_dk    = 100*(`mb_dk'-`b1_true')/`b1_true'
        local bias_pcse  = 100*(`mb_pcse'-`b1_true')/`b1_true'
    }
    else {
        local bias_naive = 100*`mb_naive'
        local bias_dk    = 100*`mb_dk'
        local bias_pcse  = 100*`mb_pcse'
    }

    post sim_res (`scen') (`Nval') (`T') (`posttrend') (1) ///
        (`n_rej_naive'/`n_ok') (`n_cov_naive'/`n_ok') (`bias_naive') (`rmse_naive') (1) ///
        (`rho1_bias_naive') (`rho1_rmse_naive') (`n_rho1_cov_naive'/`n_ok') (1)
    post sim_res (`scen') (`Nval') (`T') (`posttrend') (2) ///
        (`n_rej_dk'/`n_ok') (`n_cov_dk'/`n_ok') (`bias_dk') (`rmse_dk') (`mse_dk'/`mse_naive') ///
        (`rho1_bias_dk') (`rho1_rmse_dk') (`n_rho1_cov_dk'/`n_ok') (`rho1_mse_dk'/`rho1_mse_naive')
    post sim_res (`scen') (`Nval') (`T') (`posttrend') (3) ///
        (`n_rej_pcse'/`n_ok') (`n_cov_pcse'/`n_ok') (`bias_pcse') (`rmse_pcse') (`mse_pcse'/`mse_naive') ///
        (`rho1_bias_pcse') (`rho1_rmse_pcse') (`n_rho1_cov_pcse'/`n_ok') (`rho1_mse_pcse'/`rho1_mse_naive')

}  // posttrend
}  // T
}  // N

postclose sim_res

use "E:\ArielStuff\Simulations\sim5.dta", clear
label define scen_lbl 1 "Mild positive" 2 "Oscillatory" 3 "High persistent"
label define meth_lbl 1 "naive (oim)" 2 "xtbetark2 vce(dk)" 3 "xtbetark2 vce(pcse)"
label values scenario scen_lbl
label values method   meth_lbl
label var scenario  "Autocorrelation scenario"
label var N_panels  "Number of panels (N)"
label var T_periods "Number of time periods (T)"
label var posttrend "True post-intervention %% trend effect (0=null)"
label var method    "Estimation method"
label var reject    "Rejection rate: TYPE I ERROR when posttrend==0, POWER when posttrend!=0 -- filter on posttrend to get either"
label var coverage  "95%% CI coverage rate (slope-change coefficient)"
label var pct_bias  "Percentage bias in slope-change coefficient"
label var rmse      "Root mean squared error (slope-change coefficient)"
label var se_ratio  "Mean SE relative to naive (oim) in the same replication"
label var rho1_bias     "Percentage bias in the AR(1) coefficient (rho1)"
label var rho1_rmse     "RMSE of the AR(1) coefficient (rho1)"
label var rho1_coverage "95%% CI coverage rate for rho1"
label var rho1_se_ratio "Mean SE for rho1 relative to naive (oim) in the same replication"
save "E:\ArielStuff\Simulations\sim5.dta", replace

di as txt _newline "AR(1) primary, scen=`scen', T-subset A complete -> E:\ArielStuff\Simulations\sim5.dta"
