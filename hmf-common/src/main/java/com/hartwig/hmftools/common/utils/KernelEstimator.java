package com.hartwig.hmftools.common.utils;

/*
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*
 *    KernelEstimator.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */
public class KernelEstimator
{
    private double[] mValues;
    private double[] mWeights;
    private int mNumValues;
    private double mSumOfWeights;
    private final double mStandardDev;
    private double mPrecision;

    private static final double MAX_ERROR = 0.01;
    private static final int DEFAULT_VALUE_COUNT = 50;

    public KernelEstimator(double precision, double bandwidth)
    {
        this(precision, bandwidth, DEFAULT_VALUE_COUNT);
    }

    public KernelEstimator(double precision, double bandwidth, int expectedValues)
    {
        mPrecision = precision;
        mStandardDev = bandwidth;
        mValues = new double[expectedValues+1];
        mWeights = new double[expectedValues+1];
        mNumValues = 0;
        mSumOfWeights = 0.0;
    }

    public void addValue(double data, double weight)
    {
        if(weight != 0.0D)
        {
            data = round(data);
            int insertIndex = findNearestValue(data);
            if(mNumValues > insertIndex && mValues[insertIndex] == data)
            {
                mWeights[insertIndex] += weight;
            }
            else
            {
                if(mNumValues < mValues.length)
                {
                    int left = mNumValues - insertIndex;
                    System.arraycopy(mValues, insertIndex, mValues, insertIndex + 1, left);
                    System.arraycopy(mWeights, insertIndex, mWeights, insertIndex + 1, left);
                    mValues[insertIndex] = data;
                    mWeights[insertIndex] = weight;
                    ++mNumValues;
                }
                else
                {
                    double[] newValues = new double[mValues.length * 2];
                    double[] newWeights = new double[mValues.length * 2];
                    int left = mNumValues - insertIndex;
                    System.arraycopy(mValues, 0, newValues, 0, insertIndex);
                    System.arraycopy(mWeights, 0, newWeights, 0, insertIndex);
                    newValues[insertIndex] = data;
                    newWeights[insertIndex] = weight;
                    System.arraycopy(mValues, insertIndex, newValues, insertIndex + 1, left);
                    System.arraycopy(mWeights, insertIndex, newWeights, insertIndex + 1, left);
                    ++mNumValues;
                    mValues = newValues;
                    mWeights = newWeights;
                }
            }

            mSumOfWeights += weight;
        }
    }

    public double getProbability(double data)
    {
        double delta;
        double sum = 0.0D;
        double currentProb;
        double zLower;
        double zUpper;
        if(mNumValues == 0)
        {
            zLower = (data - mPrecision / 2.0D) / mStandardDev;
            zUpper = (data + mPrecision / 2.0D) / mStandardDev;
            return normalProbability(zUpper) - normalProbability(zLower);
        }
        else
        {
            double weightSum = 0.0D;
            int start = findNearestValue(data);

            int i;
            for(i = start; i < mNumValues; ++i)
            {
                delta = mValues[i] - data;
                zLower = (delta - mPrecision / 2.0D) / mStandardDev;
                zUpper = (delta + mPrecision / 2.0D) / mStandardDev;
                currentProb = normalProbability(zUpper) - normalProbability(zLower);
                sum += currentProb * mWeights[i];
                weightSum += mWeights[i];
                if(currentProb * (mSumOfWeights - weightSum) < sum * MAX_ERROR)
                {
                    break;
                }
            }

            for(i = start - 1; i >= 0; --i)
            {
                delta = mValues[i] - data;
                zLower = (delta - mPrecision / 2.0D) / mStandardDev;
                zUpper = (delta + mPrecision / 2.0D) / mStandardDev;
                currentProb = normalProbability(zUpper) - normalProbability(zLower);
                sum += currentProb * mWeights[i];
                weightSum += mWeights[i];
                if(currentProb * (mSumOfWeights - weightSum) < sum * MAX_ERROR)
                {
                    break;
                }
            }

            return sum / mSumOfWeights;
        }
    }

    private int findNearestValue(double key)
    {
        int low = 0;
        int high = mNumValues;

        while(low < high)
        {
            int middle = (low + high) / 2;
            double current = mValues[middle];
            if(current == key)
            {
                return middle;
            }

            if(current > key)
            {
                high = middle;
            }
            else if(current < key)
            {
                low = middle + 1;
            }
        }

        return low;
    }

    private double round(double data)
    {
        return Math.rint(data / mPrecision) * mPrecision;
    }

    // statistical methods used by the KDE
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
