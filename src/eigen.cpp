#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <cmath>
#include <float.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp ;


double line_searcher_internal (double& ss,
                               arma::mat& Co_x, // Nco * (1+p)
                               arma::colvec& Tr_total, // Ntr * 1
                               arma::colvec& coefs, // (1+p) * 1
                               arma::colvec& Newton, // (1+p) * 1
                               arma::colvec& Base_weight,
                               arma::colvec& alpha) {  
    arma::colvec weights_temp ;
    arma::mat weights_ebal ;
    arma::mat Co_x_agg ; 
    arma::colvec dif ;
    double maxdif ;
    arma::colvec penalty ;

    weights_temp = Co_x * (coefs - ss * Newton) ; // column vector
    //int n = weights_temp.n_rows ;
    //for (int i = 0; i < n; i++) {
    //  weights_temp(i, 0) = exp(weights_temp(i, 0)) ;
    //}
    weights_temp = arma::exp(weights_temp) ;
    // point wise product
    weights_ebal = weights_temp % Base_weight ;
    weights_ebal = weights_ebal.t() ;
    penalty = 2 * alpha % coefs ;

    Co_x_agg = weights_ebal * Co_x ;
    Co_x_agg = Co_x_agg.t() ;
    Co_x_agg += penalty ; 
    dif = abs(Co_x_agg - Tr_total) ;
    maxdif = dif.max() ;
    if (!isfinite(maxdif)){
        maxdif = DBL_MAX;
    }
    return(maxdif) ;

}

/*  sub function  */
// [[Rcpp::export]]
double line_searcher (arma::mat Co_x, // Nco * (1+p)
                      arma::colvec Tr_total, // Ntr * 1
                      arma::colvec coefs, // (1+p) * 1
                      arma::colvec Newton, // (1+p) * 1
                      arma::colvec Base_weight, // Nco * 1
                      arma::colvec alpha, 
                      double ss) {

    arma::mat weights_temp ;
    arma::mat weights_ebal ;
    arma::mat Co_x_agg ; 
    arma::colvec dif ; 
    double maxdif ;
    arma::colvec penalty ;

    weights_temp = Co_x * (coefs - ss * Newton) ; // column vector
    //int n = weights_temp.n_rows ;
    //for (int i = 0; i < n; i++) {
    //  weights_temp(i, 0) = exp(weights_temp(i, 0)) ;
    //}
    weights_temp = arma::exp(weights_temp) ;
    // point wise product
    weights_ebal = weights_temp % Base_weight ;
    weights_ebal = weights_ebal.t() ;
    penalty = 2 * alpha % coefs ;

    Co_x_agg = weights_ebal * Co_x ;
    Co_x_agg = Co_x_agg.t() ;
    Co_x_agg += penalty ; 
    dif = abs(Co_x_agg - Tr_total) ;
    maxdif = dif.max() ;
    if (!isfinite(maxdif)){
        maxdif = DBL_MAX;
    }
    return(maxdif) ;

}

/* optimize  line_searacher */

double Brent_fmin(double ax, 
                  double bx, 
                  double (*f)(double&, arma::mat&, arma::colvec&, arma::colvec&, arma::colvec&, arma::colvec&, arma::colvec&), 
                  arma::mat co, 
                  arma::colvec tr, 
                  arma::colvec coef, 
                  arma::colvec newton, 
                  arma::colvec base,
                  arma::colvec alpha, 
                  double tol){

    /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;

    /* Local variables */
    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;

/*  eps is approximately the square root of the relative machine precision. */
    eps = DBL_EPSILON;
    tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);

    a = ax;
    b = bx;
    v = a + c * (b - a);
    w = v;
    x = v;

    d = 0.;/* -Wall */
    e = 0.;
    fx = (*f)(x, co, tr, coef, newton, base, alpha);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;

/*  main loop starts here ----------------------------------- */

    for(;;) {
    xm = (a + b) * .5;
    tol1 = eps * fabs(x) + tol3;
    t2 = tol1 * 2.;

    /* check stopping criterion */

    if (fabs(x - xm) <= t2 - (b - a) * .5) break;
    p = 0.;
    q = 0.;
    r = 0.;
    if (fabs(e) > tol1) { /* fit parabola */

        r = (x - w) * (fx - fv);
        q = (x - v) * (fx - fw);
        p = (x - v) * q - (x - w) * r;
        q = (q - r) * 2.;
        if (q > 0.) p = -p; else q = -q;
        r = e;
        e = d;
    }

    if (fabs(p) >= fabs(q * .5 * r) ||
        p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */

        if (x < xm) e = b - x; else e = a - x;
        d = c * e;
    }
    else { /* a parabolic-interpolation step */

        d = p / q;
        u = x + d;

        /* f must not be evaluated too close to ax or bx */

        if (u - a < t2 || b - u < t2) {
        d = tol1;
        if (x >= xm) d = -d;
        }
    }

    /* f must not be evaluated too close to x */

    if (fabs(d) >= tol1)
        u = x + d;
    else if (d > 0.)
        u = x + tol1;
    else
        u = x - tol1;

    fu = (*f)(u, co, tr, coef, newton, base, alpha);

    /*  update  a, b, v, w, and x */

    if (fu <= fx) {
        if (u < x) b = x; else a = x;
        v = w;    w = x;   x = u;
        fv = fw; fw = fx; fx = fu;
    } else {
        if (u < x) a = u; else b = u;
        if (fu <= fw || w == x) {
        v = w; fv = fw;
        w = u; fw = fu;
        } else if (fu <= fv || v == x || v == w) {
        v = u; fv = fu;
        }
    }
    }
    /* end of main loop */

    return x;
}


/*  main function  */
// [[Rcpp::export]]
List hb (arma::colvec tr_total, // Ntr * 1
         arma::mat co_x, // Nco * (P+1) 
         arma::colvec coefs, // (P+1) * 1
         arma::colvec base_weight,
         arma::colvec alpha,
         int max_iterations,
         double constraint_tolerance,
         int print_level) {

    int converged = 0 ;
    int iter = 0 ;

    arma::colvec weights_temp ;
    arma::mat weights_ebal ;
    arma::mat co_x_agg ;
    arma::colvec gradient ;
    arma::colvec dif ;
    arma::mat hessian ;
    arma::colvec Coefs ;
    arma::colvec newton ;
    arma::colvec penalty ;
    arma::mat penalty_hessian ;
    arma::colvec diag;
    arma::colvec last_coefs;

    arma::colvec best_coef ; // store coef with lowest loss
    double best_maxdif ; 
    double counter = 0 ;
    double loss_new ;
    double loss_old ;
    double tol = pow(DBL_EPSILON, 0.25) ;

    double minimum ;
    double maxx = 1;
    double minn = 0.001;
//  double objective ;

    List ss_out ;


    while (iter <= max_iterations && converged == 0) {
        
        weights_temp = arma::exp(co_x * coefs) ;
        weights_ebal = weights_temp %  base_weight ; // Nco * 1
        weights_ebal = weights_ebal.t() ; // 1 * Nco

        co_x_agg = weights_ebal * co_x ; // 1 * (p+1)
        penalty = 2 * alpha % coefs;
        gradient = co_x_agg.t() - tr_total ;
        gradient += penalty ;
        dif = arma::abs(gradient) ;

        if (dif.max() < constraint_tolerance) {
            converged = 1 ;
        };


        if(isinf(dif.max())){
            coefs = last_coefs;
            weights_temp = arma::exp(co_x * coefs) ;
            weights_ebal = weights_temp %  base_weight ; // Nco * 1
            weights_ebal = weights_ebal.t() ; // 1 * Nco
            break;
        } 

        if(print_level>=2){
            Rcpp::Rcout.precision(4);
            Rcpp::Rcout<< "Iteration " << iter << ":" << " Max difference = " << dif.max() << std::fixed << std::endl;};

        diag = 2*alpha ;
        penalty_hessian = diagmat(diag) ;
        hessian = co_x.t() * (co_x.each_col() % weights_ebal.t()) + penalty_hessian ;
        last_coefs = coefs ; 

        Coefs = coefs ;
        newton = arma::solve(hessian, gradient) ;
        coefs = coefs - maxx * newton ;

        loss_new = line_searcher(co_x, tr_total, coefs, newton, base_weight, alpha, 0.0) ;
        loss_old = line_searcher(co_x, tr_total, Coefs, newton, base_weight, alpha, 0.0) ;

//        Rcout << "loss_old= " << loss_old << std::endl;

        if(print_level>=3){
            Rcout << "new loss= " << loss_new << ","<< " old loss= " << loss_old << std::endl;}

        if (loss_old <= loss_new) {
            minimum = Brent_fmin(0.0001, 1.0, &line_searcher_internal, co_x, tr_total, Coefs, newton, base_weight, alpha, tol) ;
            if(print_level>=3){Rcpp::Rcout << "LS Step Length is " << minimum << std::endl;};
            coefs = Coefs - minimum * newton ;
        }

        if (loss_old <= loss_new) {
            if (maxx <= minn){
                minn = maxx/10 ; 
            }
            minimum = Brent_fmin(minn, maxx, &line_searcher_internal, co_x, tr_total, Coefs, newton, base_weight, alpha, tol) ;
            maxx /= 2 ;
            if(print_level>=3){Rcpp::Rcout << "LS Step Length is " << minimum << std::endl;};
            coefs = Coefs - minimum * newton ;
        }

        iter = iter + 1 ;
    }

    if(converged != 1){
        coefs = best_coef ;
        weights_temp = arma::exp(co_x * coefs) ;
        weights_ebal = weights_temp %  base_weight ;
        weights_ebal = weights_ebal.t() ;
    } 

    List out ;

    out["maxdiff"] = dif.max() ;
    out["coefs"] = coefs ;
    out["Weights_ebal"] = weights_ebal.t() ;
    out["converged"] = converged ;

    return(out) ;

}