package com.hartwig.hmftools.common.utils.kde;

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
            return Statistics.normalProbability(zUpper) - Statistics.normalProbability(zLower);
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
                currentProb = Statistics.normalProbability(zUpper) - Statistics.normalProbability(zLower);
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
                currentProb = Statistics.normalProbability(zUpper) - Statistics.normalProbability(zLower);
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
}
