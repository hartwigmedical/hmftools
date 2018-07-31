package com.hartwig.hmftools.data_analyser.calcs;

import static java.lang.Math.abs;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.data_analyser.types.NmfMatrix;

public class DataUtils {

    public static double DBL_EPSILON = 1e-10;

    public static double[][] convertArray(List<List<Double>> dataSet, boolean transpose)
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

    public static boolean doublesEqual(double val1, double val2)
    {
        return abs(val1-val2) < DBL_EPSILON;
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

    public static NmfMatrix createMatrixFromListData(List<List<Double>> dataSet)
    {
        if(dataSet.isEmpty())
            return null;

        int rows = dataSet.size();
        int cols = dataSet.get(0).size();

        // create and populate the matrix
        NmfMatrix matrix = new NmfMatrix(rows, cols);

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

    public static void initRandom(NmfMatrix matrix, double min, double max, final Random rnGenerator)
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

    public void calcLeastSquares(final NmfMatrix matrix, final double[] data)
    {


    }

    public static BufferedWriter getNewFile(final String outputDir, final String fileName) throws IOException
    {
        if(outputDir.isEmpty())
            return null;

        String outputFileName = outputDir;
        if (!outputFileName.endsWith("/"))
        {
            outputFileName += "/";
        }

        outputFileName += fileName;

        Path outputFile = Paths.get(outputFileName);

        return Files.newBufferedWriter(outputFile);
    }

    public static void writeMatrixData(BufferedWriter writer, final NmfMatrix matrix, boolean asInt) throws IOException
    {
        final double[][] sigData = matrix.getData();

        for(int i = 0; i < matrix.Rows; ++i)
        {
            for(int j = 0; j < matrix.Cols; ++j) {

                if(asInt)
                    writer.write(String.format("%.0f", sigData[i][j]));
                else
                    writer.write(String.format("%.6f", sigData[i][j]));

                if(j < matrix.Cols-1)
                    writer.write(String.format(",", sigData[i][j]));
            }

            writer.newLine();
        }
    }


}
