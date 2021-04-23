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
    private double[] mValues = new double[50];
    private double[] mWeights = new double[50];
    private int mNumValues = 0;
    private double mSumOfWeights = 0.0D;
    private final double mStandardDev;
    private double mPrecision;

    private static final double MAX_ERROR = 0.01D;

    public KernelEstimator(double precision, double bandwidth)
    {
        mPrecision = precision;
        mStandardDev = bandwidth;
    }

    public void addValue(double data, double weight)
    {
        if(weight != 0.0D)
        {
            data = this.round(data);
            int insertIndex = this.findNearestValue(data);
            if(this.mNumValues > insertIndex && this.mValues[insertIndex] == data)
            {
                this.mWeights[insertIndex] += weight;
            }
            else
            {
                if(this.mNumValues < this.mValues.length)
                {
                    int left = this.mNumValues - insertIndex;
                    System.arraycopy(this.mValues, insertIndex, this.mValues, insertIndex + 1, left);
                    System.arraycopy(this.mWeights, insertIndex, this.mWeights, insertIndex + 1, left);
                    this.mValues[insertIndex] = data;
                    this.mWeights[insertIndex] = weight;
                    ++this.mNumValues;
                }
                else
                {
                    double[] newValues = new double[this.mValues.length * 2];
                    double[] newWeights = new double[this.mValues.length * 2];
                    int left = this.mNumValues - insertIndex;
                    System.arraycopy(this.mValues, 0, newValues, 0, insertIndex);
                    System.arraycopy(this.mWeights, 0, newWeights, 0, insertIndex);
                    newValues[insertIndex] = data;
                    newWeights[insertIndex] = weight;
                    System.arraycopy(this.mValues, insertIndex, newValues, insertIndex + 1, left);
                    System.arraycopy(this.mWeights, insertIndex, newWeights, insertIndex + 1, left);
                    ++this.mNumValues;
                    this.mValues = newValues;
                    this.mWeights = newWeights;
                }

            }

            this.mSumOfWeights += weight;
        }
    }

    public double getProbability(double data)
    {
        double delta;
        double sum = 0.0D;
        double currentProb;
        double zLower;
        double zUpper;
        if(this.mNumValues == 0)
        {
            zLower = (data - this.mPrecision / 2.0D) / this.mStandardDev;
            zUpper = (data + this.mPrecision / 2.0D) / this.mStandardDev;
            return Statistics.normalProbability(zUpper) - Statistics.normalProbability(zLower);
        }
        else
        {
            double weightSum = 0.0D;
            int start = this.findNearestValue(data);

            int i;
            for(i = start; i < this.mNumValues; ++i)
            {
                delta = this.mValues[i] - data;
                zLower = (delta - this.mPrecision / 2.0D) / this.mStandardDev;
                zUpper = (delta + this.mPrecision / 2.0D) / this.mStandardDev;
                currentProb = Statistics.normalProbability(zUpper) - Statistics.normalProbability(zLower);
                sum += currentProb * this.mWeights[i];
                weightSum += this.mWeights[i];
                if(currentProb * (this.mSumOfWeights - weightSum) < sum * MAX_ERROR)
                {
                    break;
                }
            }

            for(i = start - 1; i >= 0; --i)
            {
                delta = this.mValues[i] - data;
                zLower = (delta - this.mPrecision / 2.0D) / this.mStandardDev;
                zUpper = (delta + this.mPrecision / 2.0D) / this.mStandardDev;
                currentProb = Statistics.normalProbability(zUpper) - Statistics.normalProbability(zLower);
                sum += currentProb * this.mWeights[i];
                weightSum += this.mWeights[i];
                if(currentProb * (this.mSumOfWeights - weightSum) < sum * MAX_ERROR)
                {
                    break;
                }
            }

            return sum / this.mSumOfWeights;
        }
    }

    private int findNearestValue(double key)
    {
        int low = 0;
        int high = this.mNumValues;

        while(low < high)
        {
            int middle = (low + high) / 2;
            double current = this.mValues[middle];
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
        return Math.rint(data / this.mPrecision) * this.mPrecision;
    }
}
