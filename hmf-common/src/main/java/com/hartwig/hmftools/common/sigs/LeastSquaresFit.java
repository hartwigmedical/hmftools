package com.hartwig.hmftools.common.sigs;

// File:    NonNegativeLeastSquares.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.NonNegativeLeastSquares

import static com.hartwig.hmftools.common.utils.VectorUtils.copyVector;
import static com.hartwig.hmftools.common.utils.VectorUtils.initVector;
import static com.hartwig.hmftools.common.utils.MatrixUtils.copy;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Class NonNegativeLeastSquares provides a method for solving a least squares
 * minimization problem with nonnegativity constraints. The <TT>solve()</TT>
 * method finds an approximate solution to the linear system of equations
 * <B>Ax</B> = <B>b</B>, such that
 * ||<B>Ax</B>&nbsp;-&nbsp;<B>b</B>||<SUP>2</SUP> is minimized, and such that
 * <B>x</B> &gt;= <B>0</B>. The inputs to and outputs from the <TT>solve()</TT>
 * method are stored in the fields of an instance of class
 * NonNegativeLeastSquares.
 * <P>
 * The Java code is a translation of the Fortran subroutine <TT>NNLS</TT> from
 * Charles L. Lawson and Richard J. Hanson, <I>Solving Least Squares
 * Problems</I> (Society for Industrial and Applied Mathematics, 1995), page
 * 161.
 *
 * https://www.cs.rit.edu/~ark/pj/lib/edu/rit/numeric/NonNegativeLeastSquares.shtml
 *
 * @author  Alan Kaminsky
 * @version 22-Apr-2005
 */

public class LeastSquaresFit
{
    private final int M; // rows
    private final int N; // columns

    // The MxN-element A matrix for the least squares problem. On input to the solve() method, a contains the matrix A.
    // On output, a has been replaced with QA, where Q is an MxM-element orthogonal matrix generated during the solve() method's execution.
    private final double[][] mFactors;
    private final double[][] a;

    // The M-element b vector for the least squares problem. On input to the solve() method, b contains the vector b.
    // On output, b has been replaced with Qb, where Q is an MxM-element orthogonal matrix generated during the solve() method's execution.
    private final double[] mCounts;
    private final double[] b;

    // The N-element x vector for the least squares problem. On output from the solve() method, x contains the solution vector x.
    private final double[] mContribs;

    // The N-element index vector. On output from the solve() method: index[0] through index[nsetp-1] contain the indexes of the
    // elements in mContribs that are in set P, the set of positive values; that is,
    // the elements that are not forced to be zero (inactive constraints). index[nsetp] through index[N-1] contain the indexes of the
    // elements in mContribs that are in set Z, the set of zero values; that is, the elements that are forced to be zero (active constraints).
    private final int[] index;

    // The number of elements in the set P; that is, the number of positive values (inactive constraints). An output of the solve() method.
    private int nsetp;

    // The squared Euclidean norm of the residual vector, ||Ax - b||2. An output of the solve() method.
    public double normsqr;

    // Working storage.
    private final double[] w;
    private final double[] zz;
    private final double[] terms;

    // Maximum number of iterations.
    private final int itmax;

    private boolean mValid;

    // Magic numbers.
    private static final double factor = 0.01;

    private static final Logger LOGGER = LogManager.getLogger(LeastSquaresFit.class);

    public LeastSquaresFit(int rows, int cols)
    {
        if (rows <= 0 || cols <= 0)
        {
            mValid = false;
            throw new IllegalArgumentException("invalid rows or cols");
        }

        M = rows;
        N = cols;

        a = new double [rows][cols];
        mFactors = new double [rows][cols];

        b = new double [rows];
        mCounts = new double [rows];
        mContribs = new double [cols];

        index = new int [cols];

        w = new double [cols];
        zz = new double [rows];
        terms = new double [2];
        itmax = 3 * cols;
    }

    public void initialise(final double[][] factors, final double[] data)
    {
        copy(factors, a);
        copy(factors, mFactors);
        copyVector(data, b);
        copyVector(data, mCounts);

        for(int i = 0; i < index.length; ++i)
        {
            index[i] = 0;
        }

        initVector(w, 0);
        initVector(zz, 0);
        initVector(terms, 0);

        mValid = true;
    }

    public boolean valid() { return mValid; }
    public final double[] getContribs() { return mContribs; }
    public final double[] getCounts() { return mCounts; }
    public final double[][] getFactors() { return mFactors; }

    public void solve()
    {
        int i, iz, j, l, izmax, jz, jj, ip, ii;
        double sm, wmax, asave, unorm, ztest, up, alpha, t, cc, ss, temp;

        int iter = 0;

        // Initialize the arrays index and mContribs index[0] through index[nsetp-1] = set P, index[nsetp] through index[N-1] = set Z.
        for (i = 0; i < N; ++ i)
        {
            mContribs[i] = 0.0;
            index[i] = i;
        }
        nsetp = 0;

        mainloop: for (;;)
        {
            // Quit if all coefficients are already in the solution, or if M
            // columns of A have been triangularized.
            if (nsetp >= N || nsetp >= M) break mainloop;

            // Compute components of the dual (negative gradient) vector W.
            for (iz = nsetp; iz < N; ++ iz)
            {
                j = index[iz];
                sm = 0.0;
                for (l = nsetp; l < M; ++ l)
                {
                    sm += a[l][j]*b[l];
                }
                w[j] = sm;
            }

            // Find a candidate j to be moved from set Z to set P.
            candidateloop: for (;;)
            {
                // Find largest positive W[j].
                wmax = 0.0;
                izmax = -1;
                for (iz = nsetp; iz < N; ++ iz)
                {
                    j = index[iz];
                    if (w[j] > wmax)
                    {
                        wmax = w[j];
                        izmax = iz;
                    }
                }

                // If wmax <= 0, terminate. This indicates satisfaction of the
                // Kuhn-Tucker conditions.
                if (wmax <= 0.0) break mainloop;
                iz = izmax;
                j = index[iz];

                // The sign of W[j] is okay for j to be moved to set P. Begin
                // the transformation and check new diagonal element to avoid
                // near linear independence.
                asave = a[nsetp][j];
                up = constructHouseholderTransform (nsetp, nsetp+1, a, j);
                unorm = 0.0;
                for (l = 0; l < nsetp; ++ l)
                {
                    unorm += sqr (a[l][j]);
                }
                unorm = Math.sqrt (unorm);
                if (diff (unorm + Math.abs(a[nsetp][j])*factor, unorm) > 0.0)
                {
                    // Column j is sufficiently independent. Copy B into ZZ,
                    // update ZZ, and solve for ztest = proposed new value for
                    // X[j].
                    System.arraycopy (b, 0, zz, 0, M);
                    applyHouseholderTransform (nsetp, nsetp+1, a, j, up, zz);
                    ztest = zz[nsetp] / a[nsetp][j];

                    // If ztest is positive, we've found our candidate.
                    if (ztest > 0.0) break candidateloop;
                }

                // Reject j as a candidate to be moved from set Z to set P.
                // Restore a[nsetp][j], set w[j] = 0, and try again.
                a[nsetp][j] = asave;
                w[j] = 0.0;
            }

            // The index j = index[iz] has been selected to be moved from set Z
            // to set P. Update B, update indexes, apply Householder
            // transformations to columns in new set Z, zero subdiagonal
            // elements in column j, set w[j] = 0.
            System.arraycopy (zz, 0, b, 0, M);

            index[iz] = index[nsetp];
            index[nsetp] = j;
            ++ nsetp;

            jj = -1;
            for (jz = nsetp; jz < N; ++ jz)
            {
                jj = index[jz];
                applyHouseholderTransform (nsetp-1, nsetp, a, j, up, a, jj);
            }

            for (l = nsetp; l < M; ++ l)
            {
                a[l][j] = 0.0;
            }

            w[j] = 0.0;

            // Solve the triangular system. Store the solution temporarily in
            // zz.
            for (l = 0; l < nsetp; ++ l)
            {
                ip = nsetp - l;
                if (l != 0)
                {
                    for (ii = 0; ii < ip; ++ ii)
                    {
                        zz[ii] -= a[ii][jj] * zz[ip];
                    }
                }
                -- ip;
                jj = index[ip];
                zz[ip] /= a[ip][jj];
            }

            // Secondary loop begins here.
            secondaryloop: for (;;)
            {
                // Increment iteration counter.
                ++ iter;
                if (iter > itmax)
                {
                    LOGGER.error("Too many iterations");
                    return;
                }

                // See if all new constrained coefficients are feasible. If not,
                // compute alpha.
                alpha = 2.0;
                for (ip = 0; ip < nsetp; ++ ip)
                {
                    l = index[ip];
                    if (zz[ip] <= 0.0)
                    {
                        t = -mContribs[l] / (zz[ip] - mContribs[l]);
                        if (alpha > t)
                        {
                            alpha = t;
                            jj = ip;
                        }
                    }
                }

                // If all new constrained coefficients are feasible then alpha
                // will still be 2. If so, exit from secondary loop to main
                // loop.
                if (alpha == 2.0) break secondaryloop;

                // Otherwise, use alpha (which will be between 0 and 1) to
                // interpolate between the old mContribs and the new zz.
                for (ip = 0; ip < nsetp; ++ ip)
                {
                    l = index[ip];
                    mContribs[l] += alpha * (zz[ip] - mContribs[l]);
                }

                // Modify A and B and the index arrays to move coefficient i
                // from set P to set Z.
                i = index[jj];
                tertiaryloop: for (;;)
                {
                    mContribs[i] = 0.0;
                    if (jj != nsetp-1)
                    {
                        ++ jj;
                        for (j = jj; j < nsetp; ++ j)
                        {
                            ii = index[j];
                            index[j-1] = ii;
                            a[j-1][ii] =
                                    computeGivensRotation(a[j-1][ii], a[j][ii], terms);
                            a[j][ii] = 0.0;
                            cc = terms[0];
                            ss = terms[1];
                            for (l = 0; l < N; ++ l)
                            {
                                if (l != ii)
                                {
                                    // Apply Givens rotation to column l of A.
                                    temp = a[j-1][l];
                                    a[j-1][l] =  cc*temp + ss*a[j][l];
                                    a[j  ][l] = -ss*temp + cc*a[j][l];
                                }
                            }
                            // Apply Givens rotation to B.
                            temp = b[j-1];
                            b[j-1] =  cc*temp + ss*b[j];
                            b[j  ] = -ss*temp + cc*b[j];
                        }
                    }
                    -- nsetp;
                    index[nsetp] = i;

                    // See if the remaining coefficients in set P are feasible. They should be because of the way alpha was determined.
                    // If any are infeasible it is due to roundoff error. Anythat are nonpositive will be set to 0 and moved from set P to set Z.
                    for (jj = 0; jj < nsetp; ++ jj)
                    {
                        i = index[jj];
                        if (mContribs[i] <= 0.0) continue tertiaryloop;
                    }
                    break tertiaryloop;
                }

                // Copy b into zz, then solve the tridiagonal system again and continue the secondary loop
                System.arraycopy (b, 0, zz, 0, M);
                for (l = 0; l < nsetp; ++ l)
                {
                    ip = nsetp - l;
                    if (l != 0)
                    {
                        for (ii = 0; ii < ip; ++ ii)
                        {
                            zz[ii] -= a[ii][jj] * zz[ip];
                        }
                    }
                    -- ip;
                    jj = index[ip];
                    zz[ip] /= a[ip][jj];
                }
            }

            // Update mContribs from zz.
            for (ip = 0; ip < nsetp; ++ ip)
            {
                i = index[ip];
                mContribs[i] = zz[ip];
            }

            // All new coefficients are positive. Continue the main loop.
        }

        // Compute the squared Euclidean norm of the final residual vector.
        normsqr = 0.0;
        for (i = nsetp; i < M; ++ i)
        {
            normsqr += sqr (b[i]);
        }
    }

    /**
     * Construct a Householder transformation. <TT>u</TT> is an
     * <I>M</I>x<I>N</I>-element matrix used as an input and an output of this
     * method.
     *
     * @param  ipivot
     *     Index of the pivot element within the pivot vector.
     * @param  i1
     *     If <TT>i1</TT> &lt; <I>M,</I> the transformation will be constructed
     *     to zero elements indexed from <TT>i1</TT> through <I>M</I>-1. If
     *     <TT>i1</TT> &gt;= <I>M,</I> an identity transformation will be
     *     constructed.
     * @param  u
     *     An <I>M</I>x<I>N</I>-element matrix. On input, column
     *     <TT>pivotcol</TT> of <TT>u</TT> contains the pivot vector. On output,
     *     column <TT>pivotcol</TT> of <TT>u</TT>, along with the return value
     *     (<TT>up</TT>), contains the Householder transformation.
     * @param  pivotcol
     *     Index of the column of <TT>u</TT> that contains the pivot vector.
     *
     * @return
     *     The quantity <TT>up</TT> which is part of the Householder
     *     transformation.
     */
    private static double constructHouseholderTransform(int ipivot, int i1, double[][] u, int pivotcol)
    {
        int M = u.length;
        int j;
        double cl, clinv, sm, up;

        cl = Math.abs (u[ipivot][pivotcol]);

        // Construct the transformation.
        for (j = i1; j < M; ++ j)
        {
            cl = Math.max (Math.abs (u[j][pivotcol]), cl);
        }
        if (cl <= 0.0)
        {
            throw new IllegalArgumentException
                    ("NonNegativeLeastSquares.constructHouseholderTransform(): Illegal pivot vector");
        }
        clinv = 1.0 / cl;
        sm = sqr (u[ipivot][pivotcol] * clinv);
        for (j = i1; j < M; ++ j)
        {
            sm += sqr (u[j][pivotcol] * clinv);
        }
        cl = cl * Math.sqrt (sm);
        if (u[ipivot][pivotcol] > 0.0) cl = -cl;
        up = u[ipivot][pivotcol] - cl;
        u[ipivot][pivotcol] = cl;
        return up;
    }

    /**
     * Apply a Householder transformation to one column of a matrix. <TT>u</TT>
     * is an <I>M</I>x<I>N</I>-element matrix used as an input of this method.
     * <TT>c</TT> is an <I>M</I>x<I>N</I>-element matrix used as an input and
     * output of this method. <TT>ipivot</TT>, <TT>i1</TT>, <TT>u</TT>, and
     * <TT>pivotcol</TT> must be the same as in a previous call of
     * <TT>constructHouseholderTransform()</TT>, and <TT>up</TT> must be the
     * value returned by that method call.
     *
     * @param  ipivot
     *     Index of the pivot element within the pivot vector.
     * @param  i1
     *     If <TT>i1</TT> &lt; <I>M,</I> the transformation will zero elements
     *     indexed from <TT>i1</TT> through <I>M</I>-1. If <TT>i1</TT> &gt;=
     *     <I>M,</I> the transformation is an identity transformation.
     * @param  u
     *     An <I>M</I>x<I>N</I>-element matrix. On input, column
     *     <TT>pivotcol</TT> of <TT>u</TT>, along with <TT>up</TT>, contains the
     *     Householder transformation. This must be the output of a previous
     *     call of <TT>constructHouseholderTransform()</TT>.
     * @param  pivotcol
     *     Index of the column of <TT>u</TT> that contains the Householder
     *     transformation.
     * @param  up
     *     The rest of the Householder transformation. This must be the return
     *     value of the same previous call of
     *     <TT>constructHouseholderTransform()</TT>.
     * @param  c
     *     An <I>M</I>x<I>N</I>-element matrix. On input, column
     *     <TT>applycol</TT> of <TT>c</TT> contains the vector to which the
     *     Householder transformation is to be applied. On output, column
     *     <TT>applycol</TT> of <TT>c</TT> contains the transformed vector.
     * @param  applycol
     *     Index of the column of <TT>c</TT> to which the Householder
     *     transformation is to be applied.
     */
    private static void applyHouseholderTransform(int ipivot, int i1, double[][] u, int pivotcol, double up, double[][] c, int applycol)
    {
        int M = u.length;
        int i;
        double cl, b, sm;

        cl = Math.abs (u[ipivot][pivotcol]);
        if (cl <= 0.0)
        {
            throw new IllegalArgumentException
                    ("NonNegativeLeastSquares.applyHouseholderTransform(): Illegal pivot vector");
        }

        b = up * u[ipivot][pivotcol];
        // b must be nonpositive here. If b = 0, return.
        if (b == 0.0)
        {
            return;
        }
        else if (b > 0.0)
        {
            throw new IllegalArgumentException
                    ("NonNegativeLeastSquares.applyHouseholderTransform(): Illegal pivot element");
        }
        b = 1.0 / b;

        sm = c[ipivot][applycol] * up;
        for (i = i1; i < M; ++ i)
        {
            sm += c[i][applycol] * u[i][pivotcol];
        }
        if (sm != 0.0)
        {
            sm = sm * b;
            c[ipivot][applycol] += sm * up;
            for (i = i1; i < M; ++ i)
            {
                c[i][applycol] += sm * u[i][pivotcol];
            }
        }
    }

    /**
     * Apply a Householder transformation to a vector. <TT>u</TT> is an
     * <I>M</I>x<I>N</I>-element matrix used as an input of this method.
     * <TT>c</TT> is an <I>M</I>-element array used as an input and output of
     * this method. <TT>ipivot</TT>, <TT>i1</TT>, <TT>u</TT>, and
     * <TT>pivotcol</TT> must be the same as in a previous call of
     * <TT>constructHouseholderTransform()</TT>, and <TT>up</TT> must be the
     * value returned by that method call.
     *
     * @param  ipivot
     *     Index of the pivot element within the pivot vector.
     * @param  i1
     *     If <TT>i1</TT> &lt; <I>M,</I> the transformation will zero elements
     *     indexed from <TT>i1</TT> through <I>M</I>-1. If <TT>i1</TT> &gt;=
     *     <I>M,</I> the transformation is an identity transformation.
     * @param  u
     *     An <I>M</I>x<I>N</I>-element matrix. On input, column
     *     <TT>pivotcol</TT> of <TT>u</TT>, along with <TT>up</TT>, contains the
     *     Householder transformation. This must be the output of a previous
     *     call of <TT>constructHouseholderTransform()</TT>.
     * @param  pivotcol
     *     Index of the column of <TT>u</TT> that contains the Householder
     *     transformation.
     * @param  up
     *     The rest of the Householder transformation. This must be the return
     *     value of the same previous call of
     *     <TT>constructHouseholderTransform()</TT>.
     * @param  c
     *     An <I>M</I>-element array. On input, <TT>c</TT> contains the vector
     *     to which the Householder transformation is to be applied. On output,
     *     <TT>c</TT> contains the transformed vector.
     */
    private static void applyHouseholderTransform(int ipivot, int i1, double[][] u, int pivotcol, double up, double[] c)
    {
        int M = u.length;
        int i;
        double cl, b, sm;

        cl = Math.abs (u[ipivot][pivotcol]);
        if (cl <= 0.0)
        {
            throw new IllegalArgumentException
                    ("NonNegativeLeastSquares.applyHouseholderTransform(): Illegal pivot vector");
        }

        b = up * u[ipivot][pivotcol];
        // b must be nonpositive here. If b = 0, return.
        if (b == 0.0)
        {
            return;
        }
        else if (b > 0.0)
        {
            throw new IllegalArgumentException
                    ("NonNegativeLeastSquares.applyHouseholderTransform(): Illegal pivot element");
        }
        b = 1.0 / b;

        sm = c[ipivot] * up;
        for (i = i1; i < M; ++ i)
        {
            sm += c[i] * u[i][pivotcol];
        }
        if (sm != 0.0)
        {
            sm = sm * b;
            c[ipivot] += sm * up;
            for (i = i1; i < M; ++ i)
            {
                c[i] += sm * u[i][pivotcol];
            }
        }
    }

    /**
     * Compute the sine and cosine terms of a Givens rotation matrix. The terms
     * <TT>c</TT> and <TT>s</TT> are returned in <TT>terms[0]</TT> and
     * <TT>terms[1]</TT>, respectively, such that:
     * <PRE>
     *     [ c  s] * [a] = [sqrt(a^2+b^2)]
     *     [-s  c]   [b]   [      0      ]
     * </PRE>
     *
     * @param  a      Input argument.
     * @param  b      Input argument.
     * @param  terms  A 2-element array. On output, <TT>terms[0]</TT> contains
     *                <TT>c</TT> and <TT>terms[1]</TT> contains <TT>s</TT>.
     *
     * @return  sqrt(<TT>a</TT><SUP>2</SUP>+<TT>b</TT><SUP>2</SUP>).
     */
    private static double computeGivensRotation(double a, double b, double[] terms)
    {
        double xr, yr;

        if (Math.abs(a) > Math.abs(b))
        {
            xr = b/a;
            yr = Math.sqrt (1.0 + sqr (xr));
            terms[0] = sign (1.0/yr, a);
            terms[1] = terms[0]*xr;
            return Math.abs(a)*yr;
        }
        else if (b != 0.0)
        {
            xr = a/b;
            yr = Math.sqrt (1.0 + sqr (xr));
            terms[1] = sign (1.0/yr, b);
            terms[0] = terms[1]*xr;
            return Math.abs(b)*yr;
        }
        else
        {
            terms[0] = 0.0;
            terms[1] = 1.0;
            return 0.0;
        }
    }

    private static double diff(double x, double y)
    {
        return x - y;
    }

    private static double sqr (double x)
    {
        return x*x;
    }

    /**
     * Returns the number whose absolute value is x and whose sign is the same
     * as that of y. x is assumed to be nonnegative.
     */
    private static double sign(double x, double y)
    {
        return y >= 0.0 ? x : -x;
    }

}
