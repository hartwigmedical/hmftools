package com.hartwig.hmftools.common.utils;

import org.junit.Test;

public class MatrixTest
{
    @Test
    public void testMatrixCompare()
    {
        Matrix matrix1 = new Matrix(5, 5);
        Matrix matrix2 = new Matrix(5, 5);

        double[][] m1Data = matrix1.getData();
        double[][] m2Data = matrix2.getData();

        for(int i = 0; i < matrix1.Rows; ++i)
        {
            for(int j = 0; j < matrix1.Cols; ++j)
            {
                m1Data[i][j] = (i+1) * (j+1);

                if(j < matrix1.Cols-1)
                    m2Data[i][j] = m1Data[i][j];
            }
        }

        // assertTrue(NmfRun.signaturesEqual(matrix1, matrix2));
    }

}
