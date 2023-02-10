#include <cmath>
#include <iostream>

#include "vendor/cmath.h"

static double xs[] = { 1.0, 1.2, 1.5, 1.6, 1.8, 2.0 };
static double fs[] = { 5.000, 6.899, 11.180, 13.133, 18.119, 25.000 };
static const int KNOTS_COUNT = 6;

static double bCoefs[KNOTS_COUNT] = { 0.0 };
static double cCoefs[KNOTS_COUNT] = { 0.0 };
static double dCoefs[KNOTS_COUNT] = { 0.0 };

static int sevalLast = 0;


double lagrange(int n, double arrPoints[], double arrFunctions[], int m, double x);

double evalFunc(double x, double fx) {
	return fx - 6 * x - 3;
}

double lgr_quanc(double x) {
	return lagrange(KNOTS_COUNT, xs, fs, KNOTS_COUNT, x);
}

double spl_quanc(double x) {
	return seval(KNOTS_COUNT, x, xs, fs, bCoefs, cCoefs, dCoefs, &sevalLast);
}

int main(void) {
	
	std::cout << "//----LAGRANGE----" << std::endl;

	for (double x = 1; x <= 2.1; x += 0.1) {
		std::cout << "x =\t " << x << " |\t f(x) =\t " << lagrange(KNOTS_COUNT, xs, fs, KNOTS_COUNT, x) << std::endl;
	}

	double const LEFT = 1.0;
	double const RIGHT = 2.0;
	double const BISECT_STARTX = (LEFT + RIGHT) / 2;
	
	const double BISECT_ERR = 1e-5;

	double a = LEFT;
	double b = RIGHT;
	double x = 0;

	int i = 1;
	double fx = 0.0;
	double res = 0.0;
	std::cout << std::endl << "***\t CALCULATE f(x) = 6x + 3 for LAGRANGE" << std::endl;
	do {
		x = (a + b) / 2.0;
		fx = lagrange(KNOTS_COUNT, xs, fs, KNOTS_COUNT, x);
		res = evalFunc(x, fx);
		if (res >= 0) {
			b = x;
		}
		else {
			a = x;
		}
		std::cout << i << ". At x =\t " << x << ",\t f(x) - 6x - 3 =\t " << res << std::endl;
		i++;
	} while (abs(res) > BISECT_ERR);
	const double LGR_RES = x;

	std::cout << std::endl <<  "//-----SPLINE-----" << std::endl;
	int spl_flag = 0;
	spline(KNOTS_COUNT, 1, 1, 1, 1, xs, fs, bCoefs, cCoefs, dCoefs, &spl_flag);
	std::cout << "Spline interpolation flag:\t " << spl_flag << std::endl;
	
	for (double x = 1; x <= 2.1; x += 0.1) {
		std::cout << "x =\t " << x << " |\t f(x) =\t " << seval(KNOTS_COUNT, x, xs, fs, bCoefs, cCoefs, dCoefs, &sevalLast) << std::endl;
	}

	a = LEFT;
	b = RIGHT;
	x = 0;

	i = 1;
	fx = 0.0;
	res = 0.0;
	std::cout << std::endl << "***\t CALCULATE f(x) = 6x + 3 for SPLINE" << std::endl;
	do {
		x = (a + b) / 2.0;
		fx = seval(KNOTS_COUNT, x, xs, fs, bCoefs, cCoefs, dCoefs, &sevalLast);
		res = evalFunc(x, fx);
		if (res >= 0) {
			b = x;
		}
		else {
			a = x;
		}
		std::cout << i << ". At x =\t " << x << ",\t f(x) - 6x - 3 =\t " << res << std::endl;
		i++;
	} while (abs(res) > BISECT_ERR);
	const double SPL_RES = x;

	std::cout << "Lagrange evaluation result:\t " << LGR_RES << std::endl;
	std::cout << "Spline evaluation result:\t " << SPL_RES << std::endl;
	std::cout << "Absolute difference:\t " << abs(LGR_RES - SPL_RES) << std::endl;

	double lgr_integral = 0;
	double lgr_errest = 0;
	int lgr_nfe = 0;
	double lgr_posn = 0;
	int lgr_quanc_flag = 0;
	quanc8(lgr_quanc, LEFT, RIGHT, 1e-5, 1e-5, &lgr_integral, &lgr_errest, &lgr_nfe, &lgr_posn, &lgr_quanc_flag);

	double spl_integral = 0;
	double spl_errest = 0;
	int spl_nfe = 0;
	double spl_posn = 0;
	int spl_quanc_flag = 0;
	quanc8(spl_quanc, LEFT, RIGHT, 1e-5, 1e-5, &spl_integral, &spl_errest, &spl_nfe, &spl_posn, &spl_quanc_flag);

	std::cout << "Lagrange polynome integral =\t " << lgr_integral << std::endl;
	std::cout << "Lagrange polynome error estimation =\t " << lgr_errest << std::endl;
	std::cout << "Spline-function integral =\t " << spl_integral << std::endl;
	std::cout << "Spline-function error estimation =\t " << spl_errest << std::endl;
	std::cout << "Absolute difference = \t " << abs(lgr_integral - spl_integral) << std::endl;
	return 0;
}