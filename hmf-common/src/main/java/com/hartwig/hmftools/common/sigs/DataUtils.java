package com.hartwig.hmftools.common.sigs;

import static java.lang.Math.abs;
import static java.lang.Math.log10;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;

import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class DataUtils {

    public static double[] convertList(final List<Double> dataSet)
    {
        if(dataSet.isEmpty())
            return null;

        double[] dataArray = new double[dataSet.size()];

        for(int i = 0; i < dataSet.size(); ++i)
        {
            dataArray[i] = dataSet.get(i);
        }

        return dataArray;
    }

    public static double[][] convertArray(final List<List<Double>> dataSet, boolean transpose)
    {
        if(dataSet.isEmpty())
            return null;

        int rowCount = dataSet.size();
        int colCount = dataSet.get(0).size();
        double[][] dataArray = !transpose ? new double[rowCount][colCount] : new double[colCount][rowCount];

        for(int i = 0; i < rowCount; ++i)
        {
            final List<Double> data = dataSet.get(i);

            for(int j = 0; j < colCount; ++j)
            {
                if(!transpose)
                    dataArray[i][j] = data.get(j);
                else
                    dataArray[j][i] = data.get(j);
            }
        }

        return dataArray;
    }

    public static double sumVector(double[] vec)
    {
        double total = 0;
        for(int i = 0; i < vec.length; ++i)
        {
            total += vec[i];
        }

        return total;
    }

    public static void sumVectors(final double[] source, double[] dest)
    {
        if(source.length != dest.length)
            return;

        for(int i = 0; i < source.length; ++i)
        {
            dest[i] += source[i];
        }
    }

    public static double[] vectorMultiply(final double[] vec1, final double[] vec2)
    {
        if(vec1.length != vec2.length)
            return null;

        double[] output = new double[vec1.length];
        for(int i = 0; i < vec1.length; ++i)
        {
            output[i] = vec1[i] * vec2[i];
        }

        return output;
    }

    public static void initVector(double[] vec, double value)
    {
        for(int i = 0; i < vec.length; ++i)
        {
            vec[i] = value;
        }
    }

    public static boolean equalVector(double[] vec1, double[] vec2)
    {
        for(int i = 0; i < vec1.length; ++i)
        {
            if(!doublesEqual(vec1[i], vec2[i]))
                return false;
        }

        return true;
    }

    public static void vectorMultiply(double[] vec, double value)
    {
        for(int i = 0; i < vec.length; ++i)
        {
            vec[i] *= value;
        }
    }

    public static void copyVector(final double[] source, double[] dest)
    {
        if(source.length != dest.length)
            return;

        for(int i = 0; i < source.length; ++i)
        {
            dest[i] = source[i];
        }
    }

    public static void copyMatrix(final double[][] source, final double[][] dest)
    {
        if(source.length != dest.length)
            return;

        int rows = source.length;
        int cols = source[0].length;

        for(int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                dest[i][j] = source[i][j];
            }
        }
    }

    public static void addVector(final double[] source, double[] dest)
    {
        if(source.length != dest.length)
            return;

        for(int i = 0; i < source.length; ++i)
        {
            dest[i] += source[i];
        }
    }

    public static void convertToPercentages(double[] counts)
    {
        double total = sumVector(counts);

        if(total <= 0)
            return;

        for(int i = 0; i < counts.length; ++i)
        {
            counts[i] /= total;
        }
    }

    public static double capValue(double value, double minValue, double maxValue)
    {
        return max(min(value, maxValue), minValue);
    }

    public static boolean doublesEqual(double val1, double val2)
    {
        return doublesEqual(val1, val2, DBL_LARGE_EPSILON);
    }

    public static boolean doublesEqual(double val1, double val2, double epsilon)
    {
        return abs(val1-val2) < epsilon;
    }

    public static double DBL_LARGE_EPSILON = 1e-4;

    public static boolean greaterOrEqual(double val1, double val2)
    {
        return doublesEqual(val1, val2) || val1 - val2 > DBL_LARGE_EPSILON;
    }

    public static boolean greaterThan(double val1, double val2)
    {
        return val1 - val2 > DBL_LARGE_EPSILON;
    }

    public static boolean lessOrEqual(double val1, double val2)
    {
        return !greaterThan(val1, val2);
    }

    public static boolean lessThan(double val1, double val2)
    {
        return !greaterOrEqual(val1, val2);
    }

    public static int getPoissonRandom(double a, final Random rnGenerator)
    {
        double limit = Math.exp(-a);

        if(limit <= 1e-50)
            return 0;

        double prod = rnGenerator.nextDouble();

        int n = 0;
        for(; prod >= limit; n++)
        {
            prod *= rnGenerator.nextDouble();
        }
        return n;
    }

    public static int getPoissonRandomLarge(int a, Random rnGenerator)
    {
        double u = rnGenerator.nextDouble();
        double aLeft = a;
        int k = 0;
        double p = 1;
        double step = 300;

        while(true)
        {
            ++k;
            p *= rnGenerator.nextDouble();

            while(p < 1 && aLeft > 0)
            {
                if(aLeft > step)
                {
                    p *= Math.exp(step);
                    aLeft -= step;
                }
                else
                {
                    p *= Math.exp(aLeft);
                    aLeft = 0;
                }
            }


            if(p <= 1)
                break;
        }

        return k - 1;
    }

    public static SigMatrix createMatrixFromListData(List<List<Double>> dataSet)
    {
        if(dataSet.isEmpty())
            return null;

        int rows = dataSet.size();
        int cols = dataSet.get(0).size();

        // create and populate the matrix
        SigMatrix matrix = new SigMatrix(rows, cols);

        double[][] matrixData = matrix.getData();

        for(int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                matrixData[i][j] = dataSet.get(i).get(j);
            }
        }

        return matrix;
    }

    public static List<Integer> getSortedVectorIndices(final double[] data, boolean ascending)
    {
        // returns a list of indices into the original vector, being the sorted data list
        List<Integer> sortedList = Lists.newArrayList();

        for(int i = 0; i < data.length; ++i)
        {
            if(i == 0) {
                sortedList.add(i);
                continue;
            }

            int j = 0;
            for(; j < sortedList.size(); ++j)
            {
                int origIndex = sortedList.get(j);

                if(ascending && data[i] < data[origIndex])
                    break;
                else if(!ascending && data[i] > data[origIndex])
                    break;
            }

            sortedList.add(j, i);
        }

        return sortedList;
    }

    public static void initRandom(SigMatrix matrix, double min, double max, final Random rnGenerator)
    {
        // uniform random in range
        double[][] data = matrix.getData();

        for(int i= 0; i < matrix.Rows; i++)
        {
            for(int j = 0; j < matrix.Cols; j++)
            {
                data[i][j] = rnGenerator.nextDouble() * (max - min) + min;
            }
        }
    }

    public static final double[] calculateFittedCounts(final SigMatrix signatures, final double[] allocations)
    {
        double[] fittedCounts = new double[signatures.Rows];

        for(int transId = 0; transId < signatures.Cols; ++transId)
        {
            double allocation = allocations[transId];

            for(int catId = 0; catId < signatures.Rows; ++catId)
            {
                fittedCounts[catId] += allocation * signatures.get(catId, transId);
            }
        }

        return fittedCounts;
    }

    public static final int RESIDUAL_TOTAL = 0;
    public static final int RESIDUAL_PERC = 1;

    public static double[] calcResiduals(final double[] transcriptCounts, final double[] fittedCounts, double totalCounts)
    {
        double residualsTotal = 0;

        for(int catId = 0; catId < transcriptCounts.length; ++catId)
        {
            residualsTotal += abs(transcriptCounts[catId] - fittedCounts[catId]);
        }

        double residualsPerc = residualsTotal / totalCounts;

        return new double[] {residualsTotal, residualsPerc};

    }

    public static double calcLinearLeastSquares(final double[] params, final double[] data)
    {
        if(data.length != params.length)
            return 0;

        // returns the best ratio applying the params to the data assuming direct ratio (ie line through origin)
        double paramTotal = 0;
        double multTotal = 0;

        for(int i = 0; i < data.length; ++i)
        {
            paramTotal += params[i] * params[i];
            multTotal += params[i] * data[i];
        }

        return paramTotal > 0 ? multTotal/paramTotal : 0;
    }

    public static double calcAbsDiffs(final double[] set1, final double[] set2)
    {
        if(set1.length != set2.length)
            return 0;

        double diffTotal = 0;
        for(int i = 0; i < set1.length; ++i)
        {
            diffTotal += abs(set1[i] - set2[i]);
        }

        return diffTotal;
    }

    public static double calcMinPositiveRatio(final double[] params, final double[] data)
    {
        if(data.length != params.length)
            return 0;

        // returns the max ratio applying the params to the data
        // where the fit values does not exceed the actual data
        double minRatio = 0;

        for(int i = 0; i < data.length; ++i)
        {
            if(data[i] == 0)
                continue;

            double ratio = data[i] / params[i];

            if(ratio < minRatio || minRatio == 0)
                minRatio = ratio;
        }

        return minRatio;
    }

    public static String sizeToStr(double size)
    {
        return sizeToStr(size, false);
    }

    public static String doubleToStr(double size)
    {
        return sizeToStr(size, true);
    }

    public static String sizeToStr(double size, boolean withPrecision)
    {
        double log = log10(abs(size));

        if(withPrecision)
        {
            if (log >= 6)
                return String.format("%.3fM", size / 1e6);
            else if (log >= 3)
                return String.format("%.3fK", size / 1e3);
            else
                return String.format("%.3f", size);
        }
        else
        {
            if (log >= 6)
                return String.format("%.1fM", size / 1e6);
            else if (log >= 3)
                return String.format("%.1fK", size / 1e3);
            else
                return String.format("%.0f", size);
        }
    }

    public static double scaleRoundRatio(double value, int roundFactor)
    {
        int logScale = (int)round(log10(abs(value)));

        double roundUnit = pow(10, logScale - roundFactor);

        return round(value/roundUnit) * roundUnit;
    }

}
