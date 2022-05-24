package com.hartwig.hmftools.common.utils;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;

import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class MatrixUtils
{
    private static final Logger LOGGER = LogManager.getLogger(MatrixUtils.class);

    public static double sumMatrix(final double[][] data)
    {
        double total = 0;

        for(int row = 0; row < data.length; ++row)
        {
            for(int col = 0; col < data[0].length; ++col)
            {
                total += data[row][col];
            }
        }

        return total;
    }

    public static void clear(final double[][] data)
    {
        initialise(data, 0);
    }

    public static void initialise(final double[][] data, final double value)
    {
        for(int row = 0; row < data.length; ++row)
        {
            for(int col = 0; col < data[0].length; ++col)
            {
                data[row][col] = value;
            }
        }
    }

    public static void copy(final double[][] source, final double[][] dest)
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

    public static Matrix createMatrixFromListData(final List<List<Double>> dataSet)
    {
        if(dataSet.isEmpty())
            return null;

        int rows = dataSet.size();
        int cols = dataSet.get(0).size();

        // create and populate the matrix
        Matrix matrix = new Matrix(rows, cols);

        double[][] matrixData = matrix.getData();

        for(int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                matrixData[i][j] = dataSet.get(i).get(j);
            }
        }

        return matrix;
    }

    public static void scalarMultiplyRateAdjusted(final Matrix matrix, final Matrix other, double rateAdjust, int adjustColLimit)
    {
        // apply the scalar multiplication, but dampen the first X columns for the ref signatures
        final double[][] data = matrix.getData();
        final double[][] otherData = other.getData();

        for(int i = 0; i < matrix.Rows; i++)
        {
            for(int j = 0; j < matrix.Cols; j++)
            {
                if(j < adjustColLimit)
                {
                    if(otherData[i][j] >= 1)
                        data[i][j] *= 1 + (otherData[i][j] - 1) * rateAdjust;
                    else
                        data[i][j] *= 1 - (1 - otherData[i][j]) * rateAdjust;
                }
                else
                {
                    data[i][j] *= otherData[i][j];
                }
            }
        }
    }

    public static void scalarDivide(final Matrix matrix, final Matrix other)
    {
        scalarDivide(matrix, other, false);
    }

    public static void scalarDivide(final Matrix matrix, final Matrix other, boolean allowZeros)
    {
        // scalar product; this *= b
        final double[][] data = matrix.getData();
        final double[][] otherData = other.getData();

        for(int i = 0; i < matrix.Rows; i++)
        {
            for(int j = 0; j < matrix.Cols; j++)
            {
                if(otherData[i][j] == 0)
                {
                    if(allowZeros)
                        continue;

                    LOGGER.error("divide by zero at i={}, j={}", i, j);
                    return;
                }

                data[i][j] /= otherData[i][j];
            }
        }
    }

    public void scalarInvert(final Matrix matrix)
    {
        final double[][] data = matrix.getData();
        for(int i = 0; i < matrix.Rows; i++) {

            for (int j = 0; j < matrix.Cols; j++) {

                data[i][j] =  1/data[i][j];
            }
        }
    }

    public static double sum(final Matrix matrix)
    {
        return sumMatrix(matrix.getData());
    }

    public static double sumDiffSq(final Matrix matrix, final Matrix other)
    {
        // distance squared
        final double[][] data = matrix.getData();
        final double[][] otherData = other.getData();
        double d = 0;

        for(int i = 0; i < matrix.Rows; i++)
        {
            for(int j = 0; j < matrix.Cols; j++)
            {
                double v = data[i][j] - otherData[i][j];
                d += v*v;
            }
        }

        return d;
    }

    public static Matrix redimension(final Matrix other, int rows, int cols)
    {
        if(other.Rows == rows && other.Cols== cols)
            return other;

        // returns a resized matrix, copying any data that can be
        Matrix matrix = new Matrix(rows, cols);

        int minRows = min(other.Rows, matrix.Rows);
        int minCols = min(other.Cols, matrix.Cols);

        final double[][] otherData = other.getData();
        final double[][] mData = matrix.getData();

        for(int i = 0; i < minRows; i++)
        {
            for (int j = 0; j < minCols; j++)
            {
                mData[i][j] = otherData[i][j];
            }
        }

        return matrix;
    }

    public static Matrix multiply(final Matrix matrix, final Matrix other)
    {
        Matrix newMatrix = new Matrix(matrix.Rows, other.Cols);
        multiply(matrix, other, newMatrix, false);
        return newMatrix;
    }

    public static void multiply(final Matrix matrix, final Matrix other, Matrix dest, boolean initialiseDest)
    {
        // matrix multiply: c[i][j] = sum_k a[i][k] * b[k][j]
        if(matrix.Cols != other.Rows)
        {
            LOGGER.error("incorrect row or column");
            return;
        }

        final double[][] otherData = other.getData();

        if(initialiseDest)
            dest.initialise(0);

        final double[][] data = matrix.getData();
        double[][] destData = dest.getData();

        for(int i = 0; i < matrix.Rows; i++)
        {
            for(int j = 0; j < other.Cols; j++)
            {
                int commonColCount = matrix.Cols; // also = other.rowCount()

                for(int c = 0; c < commonColCount; c++)
                {
                    destData[i][j] += data[i][c] * otherData[c][j];
                }
            }
        }
    }

    public static void scalarMultiply(final Matrix matrix, final Matrix other)
    {
        // scalar product; this *= b
        final double[][] data = matrix.getData();
        final double[][] otherData = other.getData();

        for(int i = 0; i < matrix.Rows; i++)
        {
            for(int j = 0; j < matrix.Cols; j++)
            {
                data[i][j] *= otherData[i][j];
            }
        }
    }

    public static Matrix getDiff(final Matrix first, final Matrix second, boolean relative)
    {
        // return a matrix with the diffs between all entries
        Matrix results = new Matrix(first.Rows, first.Cols);

        if(first.Rows != second.Rows || first.Cols != second.Cols)
            return results;

        final double[][] fData = first.getData();
        final double[][] sData = second.getData();
        double[][] resultData = results.getData();

        for(int i = 0; i < first.Rows; i++) {

            for (int j = 0; j < first.Cols; j++) {

                double diff = fData[i][j] - sData[i][j];

                if(relative)
                {
                    if(sData[i][j] != 0)
                    {
                        resultData[i][j] = diff / sData[i][j];
                    }
                }
                else
                {
                    resultData[i][j] = diff;
                }
            }
        }

        return results;
    }

    public static Matrix extractNonZeros(final Matrix matrix)
    {
        // check for columns with all zeros and remove them from the matrix
        int nonZeroColCount = 0;

        for (int i = 0; i < matrix.Cols; ++i) {

            double colSum = sumVector(matrix.getCol(i));

            if (colSum > 0)
                ++nonZeroColCount;
        }

        if (nonZeroColCount == matrix.Cols)
            return matrix;

        Matrix newMatrix = new Matrix(matrix.Rows, nonZeroColCount);

        final double[][] mData = matrix.getData();
        double[][] nData = newMatrix.getData();

        int colIndex = 0;
        for (int i = 0; i < matrix.Cols; ++i)
        {
            double colSum = sumVector(matrix.getCol(i));

            if (colSum == 0)
                continue;

            for (int j = 0; j < matrix.Rows; j++)
            {
                nData[j][colIndex] = mData[j][i];
            }

            ++colIndex;
        }

        return newMatrix;
    }


}
