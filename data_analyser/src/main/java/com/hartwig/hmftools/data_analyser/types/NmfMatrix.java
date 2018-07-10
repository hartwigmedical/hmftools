package com.hartwig.hmftools.data_analyser.types;

import static java.lang.Math.abs;

public class NmfMatrix {

    final public int Rows;
    final public int Cols;

    private double[][] mData;

    public NmfMatrix(int r, int c)
    {
        Rows = r;
        Cols = c;

        mData = new double[r][c];
    }

    public NmfMatrix(final NmfMatrix other)
    {
        Rows = other.Rows;
        Cols = other.Cols;
        mData = new double[Rows][Cols];

        final double[][] otherData = other.getData();

        setData(otherData);
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
        for(int i = 0; i < Rows; ++i) {

            for (int j = 0; j < Cols; ++j) {

                mData[i][j] = value;
            }
        }
    }

    public int getDataCount() { return Rows * Cols; }

    public void set(int row, int col, double value)
    {
        if(row >= Rows || col >= Cols)
            return;

        mData[row][col] = value;
    }

    public double get(int row, int col)
    {
        if(row >= Rows || col >= Cols)
            return 0;

        return mData[row][col];
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
            return;

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
        // scalar product; this *= b
        final double[][] otherData = other.getData();

        for(int i = 0; i < Rows; i++)
        {
            for(int j = 0; j < Cols; j++)
            {
                if(otherData[i][j] == 0)
                    return;

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

    public double dist(NmfMatrix b)
    {
        // euclidean distance
        return Math.sqrt(sumDiffSq(b));
    }

    public boolean equals(final NmfMatrix other)
    {
        final double[][] otherData = other.getData();
        double epsilon = 0.000001;

        for(int i = 0; i < Rows; i++)
        {
            for(int j = 0; j < Cols; j++)
            {
                if(abs(mData[i][j] - otherData[i][j]) > epsilon)
                    return false;
            }
        }

        return true;
    }
}
