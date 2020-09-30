package com.hartwig.hmftools.common.stats;

public class FisherExactTest
{
    private double[] mCalcs;
    private int mCountMax;

    public FisherExactTest()
    {
        mCalcs = null;
        mCountMax = 0;
    }

    public void initialise(int maxSize)
    {
        mCalcs = new double[maxSize+1];
        mCountMax = maxSize;

        mCalcs[0] = 0.0;

        for (int i = 1; i <= mCountMax; i++)
        {
            mCalcs[i] = mCalcs[i - 1] + Math.log(i);
        }
    }

    public double calc(int withAwithB, int withANoB, int noAWithB, int noAnoB, double expectedCount)
    {
        if(withAwithB > expectedCount)
        {
            return getRightTailedP(withAwithB, noAWithB, withANoB, noAnoB);
        }
        else
        {
            return getLeftTailedP(withAwithB, noAWithB, withANoB, noAnoB);
        }
    }


    // a = with A, with B
    // b = no A, with B
    // c = with A, no B
    // d = no A, no B
    public final double getRightTailedP(int a, int b, int c, int d)
    {
        // aka 'greater than' test
        if (a + b + c + d > mCountMax)
            return Double.NaN;

        double p = 0;
        p += getP(a, b, c, d);

        int min = Math.min(c, b);
        for (int i = 0; i < min; i++)
        {
            p += getP(++a, --b, --c, ++d);
        }
        return p;
    }

    public final double getLeftTailedP(int a, int b, int c, int d)
    {
        // aka 'less than' test
        if (a + b + c + d > mCountMax)
            return Double.NaN;

        double p = 0;
        p += getP(a, b, c, d);

        int min = (a < d) ? a : d;

        for (int i = 0; i < min; i++)
        {
            double pTemp = getP(--a, ++b, ++c, --d);

            p += pTemp;
        }
        return p;
    }

    private double getP(int a, int b, int c, int d)
    {
        int n = a + b + c + d;

        if (n > mCountMax)
            return Double.NaN;

        double p = (mCalcs[a + b] + mCalcs[c + d] + mCalcs[a + c] + mCalcs[b + d]) - (mCalcs[a] + mCalcs[b] + mCalcs[c] + mCalcs[d] + mCalcs[n]);

        return Math.exp(p);
    }

    // returns small probability from right or left side
    public final double getCumulativeP(int a, int b, int c, int d)
    {
        if((a * d) >= (b * c))
        {
            return getRightTailedP(a, b, c, d);
        }
        else
        {
            return getLeftTailedP(a, b, c, d);
        }
    }

}
