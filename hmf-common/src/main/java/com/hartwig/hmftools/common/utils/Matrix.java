package com.hartwig.hmftools.common.utils;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sigs.DataUtils.doublesEqual;
import static com.hartwig.hmftools.common.utils.MatrixUtils.copy;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.utils.MatrixUtils.sumMatrix;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class Matrix
{
    final public int Rows;
    final public int Cols;

    private double[][] mData;
    private double[][] mDataTrans;

    private static final Logger LOGGER = LogManager.getLogger(Matrix.class);

    public Matrix(int r, int c)
    {
        Rows = r;
        Cols = c;

        mData = new double[r][c];
        mDataTrans = null;
    }

    public Matrix(final Matrix other)
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

        for(int i = 0; i < Rows; ++i)
        {
            for (int j = 0; j < Cols; ++j)
            {
                mDataTrans[j][i] = mData[i][j];
            }
        }
    }

    public double[][] getData() { return mData; }

    public void setData(final double[][] otherData)
    {
        copy(otherData, mData);
    }

    public void initialise(double value)
    {
        MatrixUtils.initialise(mData, value);
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

    public void setRow(int rowIndex, final int[] data)
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

    public void setCol(int colIndex, final int[] data)
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

    public Matrix transpose()
    {
        Matrix newMatrix = new Matrix(Cols, Rows);

        final double[][] otherData = newMatrix.getData();

        for(int i = 0; i < Rows; ++i) {

            for (int j = 0; j < Cols; ++j) {

                otherData[j][i] = mData[i][j];
            }
        }

        return newMatrix;
    }

    public boolean equals(final Matrix other)
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
}
