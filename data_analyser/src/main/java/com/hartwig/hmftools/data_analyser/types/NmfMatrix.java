package com.hartwig.hmftools.data_analyser.types;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.doublesEqual;

import com.hartwig.hmftools.data_analyser.calcs.DataUtils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class NmfMatrix {

    final public int Rows;
    final public int Cols;

    private double[][] mData;
    private double[][] mDataTrans;

    private static final Logger LOGGER = LogManager.getLogger(NmfMatrix.class);

    public NmfMatrix(int r, int c)
    {
        Rows = r;
        Cols = c;

        mData = new double[r][c];
        mDataTrans = null;
    }

    public NmfMatrix(final NmfMatrix other)
    {
        Rows = other.Rows;
        Cols = other.Cols;
        mData = new double[Rows][Cols];

        final double[][] otherData = other.getData();

        setData(otherData);
    }

    public void cacheTranspose()
    {
        if(mDataTrans == null)
            mDataTrans = new double[Cols][Rows];

        for(int i = 0; i < Rows; ++i) {

            for (int j = 0; j < Cols; ++j) {

                mDataTrans[j][i] = mData[i][j];
            }
        }
    }

    public double[][] getData() { return mData; }

    public void setData(final double[][] otherData)
    {
        for(int i = 0; i < Rows; ++i) {

            for (int j = 0; j < Cols; ++j) {

                mData[i][j] = otherData[i][j];
            }
        }
    }

    public void initialise(double value)
    {
        for(int i = 0; i < Rows; ++i)
        {
            for (int j = 0; j < Cols; ++j)
            {
                mData[i][j] = value;
            }
        }
    }

    public int getDataCount() { return Rows * Cols; }

    public void set(int row, int col, double value)
    {
        if(row >= Rows || col >= Cols)
        {
            LOGGER.error("incorrect row or column");
            return;
        }

        mData[row][col] = value;
    }

    public void setRow(int rowIndex, final double[] data)
    {
        if(rowIndex >= Rows)
            return;

        for(int i = 0; i < Cols; ++i)
        {
            mData[rowIndex][i] = data[i];
        }
    }

    public void setCol(int colIndex, final double[] data)
    {
        if(colIndex >= Cols)
            return;

        for(int i = 0; i < Rows; ++i)
        {
            mData[i][colIndex] = data[i];
        }
    }

    public double get(int row, int col)
    {
        if(row >= Rows || col >= Cols)
        {
            LOGGER.error("incorrect row or column");
            return 0;
        }

        return mData[row][col];
    }

    public boolean hasValidData(boolean allowNegative)
    {
        for(int i = 0; i < Rows; ++i) {

            for (int j = 0; j < Cols; ++j) {

                if(Double.isNaN(mData[i][j]) || Double.isInfinite(mData[i][j]))
                    return false;

                if(!allowNegative && mData[i][j] < 0)
                    return false;
            }
        }

        return true;
    }

    public NmfMatrix transpose()
    {
        NmfMatrix newMatrix = new NmfMatrix(Cols, Rows);

        final double[][] otherData = newMatrix.getData();

        for(int i = 0; i < Rows; ++i) {

            for (int j = 0; j < Cols; ++j) {

                otherData[j][i] = mData[i][j];
            }
        }

        return newMatrix;
    }

    public NmfMatrix multiply(final NmfMatrix other)
    {
        NmfMatrix newMatrix = new NmfMatrix(Rows, other.Cols);
        multiply(other, newMatrix, false);
        return newMatrix;
    }

    public void multiply(final NmfMatrix other, NmfMatrix dest, boolean initialiseDest)
    {
        // matrix multiply: c[i][j] = sum_k a[i][k] * b[k][j]
        if(Cols != other.Rows)
        {
            LOGGER.error("incorrect row or column");
            return;
        }

        final double[][] otherData = other.getData();

        if(initialiseDest)
            dest.initialise(0);

        double[][] destData = dest.getData();

        for(int i = 0; i < Rows; i++)
        {
            for(int j = 0; j < other.Cols; j++)
            {
                int commonColCount = Cols; // also = other.rowCount()

                for(int c = 0; c < commonColCount; c++)
                {
                    destData[i][j] += mData[i][c] * otherData[c][j];
                }
            }
        }
    }

    public void scalarMultiply(NmfMatrix other)
    {
        // scalar product; this *= b
        final double[][] otherData = other.getData();

        for(int i = 0; i < Rows; i++)
        {
            for(int j = 0; j < Cols; j++)
            {
                mData[i][j] *= otherData[i][j];
            }
        }
    }

    public void scalarMultiplyRateAdjusted(NmfMatrix other, double rateAdjust, int adjustColLimit)
    {
        // apply the scalar multiplication, but dampen the first X columns for the ref signatures
        final double[][] otherData = other.getData();

        for(int i = 0; i < Rows; i++)
        {
            for(int j = 0; j < Cols; j++)
            {
                if(j < adjustColLimit)
                {
                    if(otherData[i][j] >= 1)
                        mData[i][j] *= 1 + (otherData[i][j] - 1) * rateAdjust;
                    else
                        mData[i][j] *= 1 - (1 - otherData[i][j]) * rateAdjust;
                }
                else
                {
                    mData[i][j] *= otherData[i][j];
                }
            }
        }
    }

    public void scalarDivide(NmfMatrix other)
    {
        scalarDivide(other, false);
    }

    public void scalarDivide(NmfMatrix other, boolean allowZeros)
    {
        // scalar product; this *= b
        final double[][] otherData = other.getData();

        for(int i = 0; i < Rows; i++)
        {
            for(int j = 0; j < Cols; j++)
            {
                if(otherData[i][j] == 0)
                {
                    if(allowZeros)
                        continue;

                    LOGGER.error("divide by zero at i={}, j={}", i, j);
                    return;
                }

                mData[i][j] /= otherData[i][j];
            }
        }
    }

    public void scalarInvert()
    {
        for(int i = 0; i < Rows; i++) {

            for (int j = 0; j < Cols; j++) {

                mData[i][j] =  1/mData[i][j];
            }
        }
    }

    public double sum()
    {
        double total = 0;
        for(int i = 0; i < Rows; i++) {

            for(int j = 0; j < Cols; j++) {

                total += mData[i][j];
            }
        }

        return total;
    }

    public double[] getRow(int r)
    {
        return mData[r];
    }

    public double[] getCol(int c)
    {
        if(mDataTrans != null)
            return mDataTrans[c];

        double[] col = new double[Rows];

        for(int i = 0; i < Rows; ++i) {
            col[i] = mData[i][c];
        }

        return col;
    }

    public double sumDiffSq(final NmfMatrix other)
    {
        // distance squared
        final double[][] otherData = other.getData();
        double d = 0;

        for(int i = 0; i < Rows; i++)
        {
            for(int j = 0; j < Cols; j++)
            {
                double v = mData[i][j] - otherData[i][j];
                d += v*v;
            }
        }

        return d;
    }

    public static NmfMatrix redimension(final NmfMatrix other, int rows, int cols)
    {
        if(other.Rows == rows && other.Cols== cols)
            return other;

        // returns a resized matrix, copying any data that can be
        NmfMatrix matrix = new NmfMatrix(rows, cols);

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

    public double euclideanDist(NmfMatrix b)
    {
        // euclidean distance
        return Math.sqrt(sumDiffSq(b));
    }

    public boolean equals(final NmfMatrix other)
    {
        final double[][] otherData = other.getData();

        for(int i = 0; i < Rows; i++)
        {
            for(int j = 0; j < Cols; j++)
            {
                if(!doublesEqual(mData[i][j], otherData[i][j]))
                    return false;
            }
        }

        return true;
    }

    public static NmfMatrix getDiff(final NmfMatrix first, final NmfMatrix second, boolean relative)
    {
        // return a matrix with the diffs between all entries
        NmfMatrix results = new NmfMatrix(first.Rows, first.Cols);

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

    public static NmfMatrix extractNonZeros(final NmfMatrix matrix) {

        // check for columns with all zeros and remove them from the matrix
        int nonZeroColCount = 0;

        for (int i = 0; i < matrix.Cols; ++i) {

            double colSum = DataUtils.sumVector(matrix.getCol(i));

            if (colSum > 0)
                ++nonZeroColCount;
        }

        if (nonZeroColCount == matrix.Cols)
            return matrix;

        NmfMatrix newMatrix = new NmfMatrix(matrix.Rows, nonZeroColCount);

        final double[][] mData = matrix.getData();
        double[][] nData = newMatrix.getData();

        int colIndex = 0;
        for (int i = 0; i < matrix.Cols; ++i) {

            double colSum = DataUtils.sumVector(matrix.getCol(i));

            if (colSum == 0)
                continue;

            for (int j = 0; j < matrix.Rows; j++) {
                nData[j][colIndex] = mData[j][i];
            }

            ++colIndex;
        }

        return newMatrix;
    }
}
