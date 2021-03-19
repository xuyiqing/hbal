#include <RcppEigen.h>
#include <cmath>
#include <float.h>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp ;

/*  sub function  */
double line_searcher_internal (double& ss,
                               Eigen::MatrixXd& Co_x, // Nco * (1+p)
	                           Eigen::VectorXd& Tr_total, // Ntr * 1
	                           Eigen::VectorXd& coefs, // (1+p) * 1
	                           Eigen::VectorXd& Newton, // (1+p) * 1
	                           Eigen::VectorXd& Base_weight, // Nco * 1
	                           Eigen::VectorXd& alpha) {  
	Eigen::VectorXd weights_temp ;
	Eigen::VectorXd weights_ebal ;
	Eigen::VectorXd Co_x_agg ; 
	Eigen::VectorXd dif ; 
	Eigen::VectorXd penalty ;
	double maxdif ;

	weights_temp = (Co_x * (coefs - ss * Newton)).array().exp().matrix() ;
    // point wise product
	weights_ebal = (weights_temp.array() * Base_weight.array()).matrix() ;
	//weights_ebal = weights_ebal.transpose() ;
	penalty = (2 * alpha.array() * coefs.array()).matrix() ;

	Co_x_agg = Co_x.transpose() * weights_ebal;
    //Co_x_agg = Co_x_agg.transpose() ;
	Co_x_agg += penalty ; 
	dif = (Co_x_agg - Tr_total).cwiseAbs() ;
	maxdif = dif.array().maxCoeff() ;
    if (std::isnan(maxdif) or std::isinf(maxdif)){
        maxdif = DBL_MAX;
    }
	return(maxdif) ;

}


/*  sub function  */
// [[Rcpp::export]]
double line_searcher (Eigen::MatrixXd Co_x, // Nco * (1+p)
	                  Eigen::VectorXd Tr_total, // Ntr * 1
	                  Eigen::VectorXd coefs, // (1+p) * 1
	                  Eigen::VectorXd Newton, // (1+p) * 1
	                  Eigen::VectorXd Base_weight, // Nco * 1
	                  Eigen::VectorXd alpha,
	                  double ss) {  
	Eigen::VectorXd weights_temp ;
	Eigen::VectorXd weights_ebal ;
	Eigen::VectorXd Co_x_agg ; 
	Eigen::VectorXd dif ; 
	Eigen::VectorXd penalty ;
	double maxdif ;

	weights_temp = (Co_x * (coefs - ss * Newton)).array().exp().matrix() ; // column vector
	
    // point wise product
    weights_ebal = (weights_temp.array() * Base_weight.array()).matrix() ;
	//weights_ebal = weights_ebal.transpose() ;
	penalty = (2 * alpha.array() * coefs.array()).matrix() ;

    Co_x_agg = Co_x.transpose() * weights_ebal;
	//Co_x_agg = Co_x_agg.transpose() ;
    Co_x_agg = Co_x_agg + penalty ; 
	dif = (Co_x_agg - Tr_total).cwiseAbs() ;
	maxdif = dif.array().maxCoeff() ;
    if (std::isnan(maxdif) or std::isinf(maxdif)){
        maxdif = DBL_MAX;
    }
	return(maxdif) ;
}

/*  sub function  */
/*  barebone implementation of R's optimize(), adapted from optimize.c
    https://github.com/SurajGupta/r-source/blob/master/src/library/stats/src/optimize.c  */

double Brent_fmin(double ax, 
                  double bx, 
                  double (*f)(double&, 
                              Eigen::MatrixXd&, 
                              Eigen::VectorXd&, 
                              Eigen::VectorXd&, 
                              Eigen::VectorXd&, 
                              Eigen::VectorXd&, 
                              Eigen::VectorXd&), 
                  Eigen::MatrixXd co, 
                  Eigen::VectorXd tr, 
                  Eigen::VectorXd coef, 
                  Eigen::VectorXd newton, 
                  Eigen::VectorXd base, 
                  Eigen::VectorXd alpha, 
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
List hb (Eigen::VectorXd tr_total, // Ntr * 1
	     Eigen::MatrixXd co_x, // Nco * (P+1) 
	     Eigen::VectorXd coefs, // (P+1) * 1
	     Eigen::VectorXd base_weight,
	     Eigen::VectorXd alpha,
	     int max_iterations,
	     double constraint_tolerance,
	     int print_level) {

	int converged = 0 ;
	int iter = 0 ;

	Eigen::VectorXd weights_temp ;
	Eigen::VectorXd weights_ebal ;
	Eigen::VectorXd co_x_agg ;
	Eigen::VectorXd gradient ;
	Eigen::VectorXd dif ;
	Eigen::MatrixXd hessian ;
	Eigen::VectorXd Coefs ;
	Eigen::VectorXd newton ;
	Eigen::VectorXd penalty ;
	//mat diag_penalty;
	Eigen::MatrixXd penalty_hessian ;
	Eigen::VectorXd diag;
	//Eigen::VectorXd w;
    Eigen::VectorXd last_coefs;
	//double counter = 0 ;
	double loss_new ;
	double loss_old ;
	double minimum ;
    double tol = pow(DBL_EPSILON, 0.25) ;
    double err_tol = 1e-6 ;
    double maxx = 1;
    double minn = 0.001;
//  double objective ;
	List ss_out ;


	while (iter <= max_iterations && converged == 0) {
        
        weights_temp = (co_x * coefs).array().exp().matrix() ;
        
        weights_ebal = (weights_temp.array() * base_weight.array()).matrix() ;
        //w = weights_ebal.transpose();
        //weights_ebal = weights_ebal.transpose() ; // 1 * Nco

        co_x_agg = co_x.transpose() * weights_ebal; // 1 * (p+1)
		penalty = (2 * alpha.array() * coefs.array()).matrix() ;
		//co_x_agg = co_x_agg.transpose() ;
        gradient = co_x_agg - tr_total ;
        gradient +=  penalty ;
        dif = (gradient).cwiseAbs() ;
        if (dif.array().maxCoeff() < constraint_tolerance) {
        	converged = 1 ;
        };

        if(std::isinf(dif.array().maxCoeff()) or std::isnan(dif.array().maxCoeff())){
            coefs = last_coefs;
            weights_temp = (co_x * coefs).array().exp().matrix() ;
            weights_ebal = (weights_temp.array() * base_weight.array()).matrix() ;
            weights_ebal = weights_ebal.transpose() ;
            break;
        } 

        if(print_level>=2){
        	Rcpp::Rcout.precision(4);
        	Rcpp::Rcout<< "Iteration " << iter << ":" << " Max difference = " << dif.array().maxCoeff() << std::fixed << std::endl;};

        diag = (2*alpha).transpose() ;
        penalty_hessian = (diag).asDiagonal() ;
        hessian = co_x.transpose() * (co_x.array().colwise() * weights_ebal.array()).matrix() + penalty_hessian ;	

        last_coefs = coefs ;      

        Coefs = coefs ;
        newton = hessian.llt().solve(gradient) ;
        double relative_error = (hessian*newton - gradient).norm() / newton.norm();
//        Rcout<< "error =  " << (hessian*newton - gradient).norm() << "; " << " relative error = " << relative_error << std::endl;
        if (relative_error > err_tol or std::isnan(relative_error)){
            newton = hessian.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(gradient) ;
//            Rcout<<"lldterror =  " << (hessian*newton - gradient).norm() << "; " << "lldt error =  " << (hessian*newton - gradient).norm()/ newton.norm()<< std::endl;
        }
        coefs = coefs - maxx * newton ;

        loss_new = line_searcher(co_x, tr_total, coefs, newton, base_weight, alpha, 0.0) ;
        loss_old = line_searcher(co_x, tr_total, Coefs, newton, base_weight, alpha, 0.0) ;

        if(std::isinf(loss_new) or std::isnan(loss_new)){
        	coefs = Coefs;
        	break;
        }

        if(print_level>=3){
        	Rcout << "new loss= " << loss_new << ","<< " old loss= " << loss_old << std::endl;}

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

	List out ;

	out["maxdiff"] = dif.array().maxCoeff() ;
	out["coefs"] = coefs ;
	out["Weights_ebal"] = weights_ebal ;
	out["converged"] = converged ;

	return(out) ;

}
