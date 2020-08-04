package com.hartwig.hmftools.common.sigs;

import static java.lang.Math.abs;
import static java.lang.Math.log10;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.sigs.VectorUtils.sumVector;

import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;

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

    public static SigMatrix createMatrixFromListData(final List<List<Double>> dataSet)
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
