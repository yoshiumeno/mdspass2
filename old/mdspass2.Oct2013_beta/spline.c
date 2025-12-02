/* spline.f -- translated by f2c (version 20090411).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/*   SPLINE FUNCTIONS */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* Subroutine */ int splc_(double *x, int *n, double *y, 
	double *df, int *iopt, double *c__, int *nc, int *ier)
{
    /* System generated locals */
    int c_dim1, c_offset, i__1, i__2;

    /* Local variables */
    static double d__[2], h__;
    static int i__, j, k;
    static double a1, a2, h1;
    static int i1, i2;
    static double h2, x1, x2, x3, y1, y2, ec[4];
    static int ib, ke, ii;
    static double hh;
    static int ks;
    static double hy, dy1, dy2, piv;
    static int ider;

/* ************************************************************************ */
/* *  COMPUTE THE COEFFICIENTS OF THE CUBIC SPLINE.                       * */
/* *  PARAMETERS                                                          * */
/* *    (1) X: 1-DIM. ARRAY FOR KNOWN POINTS                              * */
/* *    (2) N: NUMBER OF KNOWN POINTS                                     * */
/* *    (3) Y: 1-DIM. ARRAY FOR FUNCTION'S VALUES ON KNOWN POINTS         * */
/* *    (4) DF: 1-DIM. ARRAY FOR DIFFERENTIALS AT END POINTS              * */
/* *    (5) IOPT: 1-DIM. ARRAY SPECIFYING THE CONTENT OF DF               * */
/* *    (6) C: 2-DIM. WORKING ARRAY                                       * */
/* *    (7) NC: ROW SIZE OF THE ARRAY (C)                                 * */
/* *    (8) IER: ERROR CODE                                               * */
/* *  COPYRIGHT   T. OGUNI   JUNE 30 1989    VERSION 1.0                  * */
/* ************************************************************************ */
/* C */
    /* Parameter adjustments */
    --y;
    --x;
    --df;
    --iopt;
    c_dim1 = *nc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */
    if (*n < 2 || *nc < *n - 1 || iopt[1] < 1 || iopt[1] > 3 || iopt[2] < 1 ||
	     iopt[2] > 3) {
	*ier = 2;
/*       WRITE(*,*) '(SUBR. SPLC) INVALID ARGUMENT.',N,NC,IOPT(1),IOPT(2) */
	return 0;
    }
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (x[i__] >= x[i__ + 1]) {
	    *ier = 1;
/*        WRITE(*,*) '(SUBR. SPLC) X SHOULD SATISFY UPWARD ORDER.' */
	    return 0;
	}
/* L5: */
    }
    *ier = 0;
/* C  SET THE END CONDITIONS. */
    ii = 2;
    ks = 1;
    ke = min(4,*n);
    ider = 1;
    for (i__ = 1; i__ <= 2; ++i__) {
	i1 = (i__ << 1) - 1;
	i2 = i__ << 1;
	ib = iopt[i__];
/*       GO TO (10, 20, 30), IB */
	if (ib == 1) {
	    goto L10;
	} else if (ib == 2) {
	    goto L20;
	} else {
	    goto L30;
	}
L10:
	ec[i1 - 1] = 0.;
	ec[i2 - 1] = df[i__] * 2.;
	goto L70;
L20:
	d__[i__ - 1] = df[i__];
L25:
	if (i__ == 2) {
	    ii = *n;
	}
	h__ = x[ii] - x[ii - 1];
	ec[i1 - 1] = 1.;
	hy = y[ii] - y[ii - 1];
	ec[i2 - 1] = (hy / h__ - d__[i__ - 1]) * 6. / h__;
	if (i__ == 2) {
	    ec[i2 - 1] = -ec[i2 - 1];
	}
	goto L70;
L30:
	if (i__ != 1) {
/* Computing MAX */
	    i__1 = 1, i__2 = *n - 3;
	    ks = max(i__1,i__2);
	    ke = *n;
	    ider = *n;
	}
	a2 = 0.;
	d__[i__ - 1] = 0.;
	i__1 = ke;
	for (k = ks; k <= i__1; ++k) {
	    if (ider != k) {
		a1 = 1.;
		i__2 = ke;
		for (j = ks; j <= i__2; ++j) {
		    if (j != ider && j != k) {
			x1 = x[ider] - x[j];
			x2 = x[k] - x[j];
			a1 = a1 * x1 / x2;
		    }
/* L50: */
		}
		x3 = x[k] - x[ider];
		d__[i__ - 1] += a1 * y[k] / x3;
		a2 -= 1. / x3;
	    }
/* L60: */
	}
	d__[i__ - 1] += y[ider] * a2;
	goto L25;
L70:
	;
    }
/* C  SET THE ELEMENTS FOR THE SYMMETRIC TRIDIAGONAL EQUATION. */
    if (*n != 2) {
	h1 = x[2] - x[1];
	y1 = y[2] - y[1];
	i__1 = *n - 1;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    h2 = x[i__ + 1] - x[i__];
	    y2 = y[i__ + 1] - y[i__];
	    hh = h1 + h2;
	    c__[i__ + c_dim1] = h2 / hh;
	    c__[i__ + (c_dim1 << 1)] = 1. - c__[i__ + c_dim1];
	    c__[i__ + c_dim1 * 3] = (y2 / h2 - y1 / h1) * 6. / hh;
	    h1 = h2;
	    y1 = y2;
	}
    }
/* C  SOLVE THE EQUATION */
    c__[c_dim1 + 1] = -ec[0] * .5;
    c__[(c_dim1 << 1) + 1] = ec[1] * .5;
    if (*n != 2) {
	i__1 = *n - 1;
	for (k = 2; k <= i__1; ++k) {
	    piv = c__[k + (c_dim1 << 1)] * c__[k - 1 + c_dim1] + 2.;
	    c__[k + c_dim1] = -c__[k + c_dim1] / piv;
	    c__[k + (c_dim1 << 1)] = (c__[k + c_dim1 * 3] - c__[k + (c_dim1 <<
		     1)] * c__[k - 1 + (c_dim1 << 1)]) / piv;
	}
    }
    dy1 = (ec[3] - ec[2] * c__[*n - 1 + (c_dim1 << 1)]) / (ec[2] * c__[*n - 1 
	    + c_dim1] + 2.);
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = *n - i__;
	dy2 = c__[k + c_dim1] * dy1 + c__[k + (c_dim1 << 1)];
	h__ = x[k + 1] - x[k];
	c__[k + c_dim1 * 3] = (dy1 - dy2) / (h__ * 6.);
	c__[k + (c_dim1 << 1)] = dy2 * .5;
	c__[k + c_dim1] = (y[k + 1] - y[k]) / h__ - (c__[k + (c_dim1 << 1)] + 
		c__[k + c_dim1 * 3] * h__) * h__;
	dy1 = dy2;
    }
/* C */
    return 0;
} /* splc_ */

/* Subroutine */ int splf_(double *x, int *n, double *y, 
	double *c__, int *nc, double *v, int *m, double *
	f, int *ier)
{
    /* System generated locals */
    int c_dim1, c_offset, i__1;

    /* Local variables */
    static int i__, k;
    static double v1, v2;

/* ************************************************************************ */
/* *  INTERPOLATION BY THE CUBIC SPLINE.                                  * */
/* *  PARAMETERS                                                          * */
/* *    (1) X: 1-DIM. ARRAY FOR KNOWN POINTS                              * */
/* *    (2) N: NUMBER OF KNOWN POINTS                                     * */
/* *    (3) Y: 1-DIM. ARRAY FOR FUNCTION'S VALUES ON KNOWN POINTS         * */
/* *    (4) C: 2-DIM. WORKING ARRAY                                       * */
/* *    (5) NC: ROW SIZE OF THE ARRAY (C)                                 * */
/* *    (6) V: 1-DIM. ARRAY FOR POINTS WHICH INTERPOLATION MUST BE MADE   * */
/* *    (7) M: NUMBER OF POINTS FOR WHICH INTERPOLATION MUST BE MADE      * */
/* *    (8) F: 1-DIM. WORKING ARRAY                                       * */
/* *    (9) IER: ERROR CODE                                               * */
/* *  COPYRIGHT   T. OGUNI   JUNE 30 1989   VERSION 1.0                   * */
/* ************************************************************************ */
/* C */
    /* Parameter adjustments */
    --y;
    --x;
    c_dim1 = *nc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --f;
    --v;

    /* Function Body */
    if (*n < 2 || *m < 1 || *nc < *n - 1) {
	*ier = 2;
/*       WRITE(*,*) '(SUBR. SPLF) INVALID ARGUMENT. ', N, NC, M */
	return 0;
    }
    *ier = 0;
    i__ = 1;
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	v1 = v[k] - x[i__];
/*       IF (V1) 10, 30, 40 */
	if (v1 < 0.f) {
	    goto L10;
	} else if (v1 == 0.f) {
	    goto L30;
	} else {
	    goto L40;
	}
L10:
	if (i__ > 1) {
	    goto L20;
	}
	*ier = 1;
	goto L80;
L20:
	--i__;
	v1 = v[k] - x[i__];
/*       IF (V1) 10, 30, 80 */
	if (v1 < 0.f) {
	    goto L10;
	} else if (v1 == 0.f) {
	    goto L30;
	} else {
	    goto L80;
	}
L30:
	f[k] = y[i__];
	goto L90;
L40:
	if (i__ < *n) {
	    goto L50;
	}
	*ier = 1;
	i__ = *n - 1;
	goto L80;
L50:
	v2 = v[k] - x[i__ + 1];
/*       IF (V2) 80, 60, 70 */
	if (v2 < 0.f) {
	    goto L80;
	} else if (v2 == 0.f) {
	    goto L60;
	} else {
	    goto L70;
	}
L60:
	++i__;
	goto L30;
L70:
	++i__;
	v1 = v2;
	goto L40;
L80:
	f[k] = y[i__] + v1 * (c__[i__ + c_dim1] + v1 * (c__[i__ + (c_dim1 << 
		1)] + v1 * c__[i__ + c_dim1 * 3]));
L90:
	;
    }
/* C */
    return 0;
} /* splf_ */

/* Subroutine */ int spld_(double *x, int *n, double *c__, 
	int *nc, double *v, int *m, double *d1, double *
	d2, int *ier)
{
    /* System generated locals */
    int c_dim1, c_offset, i__1;

    /* Local variables */
    static int i__, k;
    static double t, v1, v2;

/* ************************************************************************ */
/* *  DIFFERENTIATION BY THE CUBIC SPLINE.                                * */
/* *  PARAMETERS                                                          * */
/* *    (1) X: 1-DIM. ARRAY FOR KNOWN POINTS                              * */
/* *    (2) N: NUMBER OF KNOWN POINTS                                     * */
/* *    (3) C: 2-DIM. WORKING ARRAY                                       * */
/* *    (4) NC: ROW SIZE OF THE ARRAY (C)                                 * */
/* *    (5) V: 1-DIM. ARRAY FOR POINTS WHICH INTERPOLATION MUST BE MADE   * */
/* *    (6) M: NUMBER OF POINTS FOR WHICH INTERPOLATION MUST BE MADE      * */
/* *    (7) D1: 1-DIM. ARRAY FOR FIRST ORDER DIFFERENTIALS                * */
/* *    (8) D2: 1-DIM. ARRAY FOR SECOND ORDER DIFFERENTIALS               * */
/* *    (9) IER: ERROR CODE                                               * */
/* *  COPYRIGHT   T. OGUNI   JUNE 30 1989    VERSION 1.0                  * */
/* ************************************************************************ */
/* C */
    /* Parameter adjustments */
    --x;
    c_dim1 = *nc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --d2;
    --d1;
    --v;

    /* Function Body */
    if (*n < 2 || *nc < *n - 1 || *m < 1) {
	*ier = 2;
/*       WRITE(*,*) '(SUBR. SPLD) INVALID ARGUMENT. ', N, NC, M */
	return 0;
    }
    *ier = 0;
    i__ = 1;
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	v1 = v[k] - x[i__];
/*       IF (V1) 10, 30, 40 */
	if (v1 < 0.f) {
	    goto L10;
	} else if (v1 == 0.f) {
	    goto L30;
	} else {
	    goto L40;
	}
L10:
	if (i__ > 1) {
	    goto L20;
	}
	*ier = 1;
	goto L80;
L20:
	--i__;
	v1 = v[k] - x[i__];
/*       IF (V1) 10, 30 ,80 */
	if (v1 < 0.f) {
	    goto L10;
	} else if (v1 == 0.f) {
	    goto L30;
	} else {
	    goto L80;
	}
L30:
	d1[k] = c__[i__ + c_dim1];
	d2[k] = c__[i__ + (c_dim1 << 1)] + c__[i__ + (c_dim1 << 1)];
	goto L90;
L40:
	if (i__ < *n) {
	    goto L50;
	}
	*ier = 1;
	i__ = *n - 1;
	goto L80;
L50:
	v2 = v[k] - x[i__ + 1];
/*       IF (V2) 80, 60, 70 */
	if (v2 < 0.f) {
	    goto L80;
	} else if (v2 == 0.f) {
	    goto L60;
	} else {
	    goto L70;
	}
L60:
	if (i__ >= *n - 1) {
	    goto L80;
	}
	++i__;
	goto L30;
L70:
	++i__;
	v1 = v2;
	goto L40;
L80:
	t = c__[i__ + (c_dim1 << 1)] + c__[i__ + c_dim1 * 3] * 3. * v1;
	d1[k] = c__[i__ + c_dim1] + (t + c__[i__ + (c_dim1 << 1)]) * v1;
	d2[k] = t + t;
L90:
	;
    }
    return 0;
} /* spld_ */

