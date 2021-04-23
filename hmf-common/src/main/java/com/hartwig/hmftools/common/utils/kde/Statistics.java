//
// Source code recreated from a .class file by IntelliJ IDEA
// (powered by Fernflower decompiler)
//

package com.hartwig.hmftools.common.utils.kde;

/**
 * Class implementing some distributions, tests, etc. The code is mostly adapted from the CERN
 * Jet Java libraries:
 * <p>
 * Copyright 2001 University of Waikato
 * Copyright 1999 CERN - European Organization for Nuclear Research.
 * Permission to use, copy, modify, distribute and sell this software and its documentation for
 * any purpose is hereby granted without fee, provided that the above copyright notice appear
 * in all copies and that both that copyright notice and this permission notice appear in
 * supporting documentation.
 * CERN and the University of Waikato make no representations about the suitability of this
 * software for any purpose. It is provided "as is" without expressed or implied warranty.
 *
 * @author peter.gedeck@pharma.Novartis.com
 * @author wolfgang.hoschek@cern.ch
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @author Richard Kirkby (rkirkby@cs.waikato.ac.nz)
 * @version $Revision: 5619 $
 */
class Statistics
{
    protected static final double[] P1 =
            new double[] { 4.0554489230596245D, 31.525109459989388D, 57.16281922464213D, 44.08050738932008D, 14.684956192885803D,
                    2.1866330685079025D, -0.1402560791713545D, -0.03504246268278482D, -8.574567851546854E-4D };

    public static double normalProbability(double a)
    {
        double x = a * 0.7071067811865476D;
        double z = Math.abs(x);
        double y;
        if(z < 0.7071067811865476D)
        {
            y = 0.5D + 0.5D * errorFunction(x);
        }
        else
        {
            y = 0.5D * errorFunctionComplemented(z);
            if(x > 0.0D)
            {
                y = 1.0D - y;
            }
        }

        return y;
    }

    public static double errorFunction(double x)
    {
        double[] T = new double[] { 9.604973739870516D, 90.02601972038427D, 2232.005345946843D, 7003.325141128051D, 55592.30130103949D };
        double[] U = new double[] { 33.56171416475031D, 521.3579497801527D, 4594.323829709801D, 22629.000061389095D, 49267.39426086359D };
        if(Math.abs(x) > 1.0D)
        {
            return 1.0D - errorFunctionComplemented(x);
        }
        else
        {
            double z = x * x;
            double y = x * polevl(z, T, 4) / p1evl(z, U, 5);
            return y;
        }
    }

    public static double errorFunctionComplemented(double a)
    {
        double[] P = new double[] { 2.461969814735305E-10D, 0.5641895648310689D, 7.463210564422699D, 48.63719709856814D, 196.5208329560771D,
                526.4451949954773D, 934.5285271719576D, 1027.5518868951572D, 557.5353353693994D };
        double[] Q = new double[] { 13.228195115474499D, 86.70721408859897D, 354.9377788878199D, 975.7085017432055D, 1823.9091668790973D,
                2246.3376081871097D, 1656.6630919416134D, 557.5353408177277D };
        double[] R = new double[] { 0.5641895835477551D, 1.275366707599781D, 5.019050422511805D, 6.160210979930536D, 7.4097426995044895D,
                2.9788666537210022D };
        double[] S = new double[] { 2.2605286322011726D, 9.396035249380015D, 12.048953980809666D, 17.08144507475659D, 9.608968090632859D,
                3.369076451000815D };
        double x;
        if(a < 0.0D)
        {
            x = -a;
        }
        else
        {
            x = a;
        }

        if(x < 1.0D)
        {
            return 1.0D - errorFunction(a);
        }
        else
        {
            double z = -a * a;
            if(z < -709.782712893384D)
            {
                return a < 0.0D ? 2.0D : 0.0D;
            }
            else
            {
                z = Math.exp(z);
                double p;
                double q;
                if(x < 8.0D)
                {
                    p = polevl(x, P, 8);
                    q = p1evl(x, Q, 8);
                }
                else
                {
                    p = polevl(x, R, 5);
                    q = p1evl(x, S, 6);
                }

                double y = z * p / q;
                if(a < 0.0D)
                {
                    y = 2.0D - y;
                }

                return y == 0.0D ? (a < 0.0D ? 2.0D : 0.0D) : y;
            }
        }
    }

    public static double p1evl(double x, double[] coef, int N)
    {
        double ans = x + coef[0];

        for(int i = 1; i < N; ++i)
        {
            ans = ans * x + coef[i];
        }

        return ans;
    }

    public static double polevl(double x, double[] coef, int N)
    {
        double ans = coef[0];

        for(int i = 1; i <= N; ++i)
        {
            ans = ans * x + coef[i];
        }

        return ans;
    }
}
