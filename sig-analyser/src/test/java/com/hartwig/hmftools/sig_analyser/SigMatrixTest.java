package com.hartwig.hmftools.sig_analyser;

import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.sig_analyser.common.SigMatrix;
import com.hartwig.hmftools.sig_analyser.nmf.NmfRun;

import org.junit.Test;

public class SigMatrixTest
{
    @Test
    public void testSigCompare()
    {
        SigMatrix sigs1 = new SigMatrix(5, 5);
        SigMatrix sigs2 = new SigMatrix(5, 5);

        double[][] s1Data = sigs1.getData();
        double[][] s2Data = sigs2.getData();

        for(int i = 0; i < sigs1.Rows; ++i)
        {
            for(int j = 0; j < sigs1.Cols; ++j)
            {
                s1Data[i][j] = (i+1) * (j+1);

                if(j < sigs1.Cols-1)
                    s2Data[i][j] = s1Data[i][j];
            }
        }

        // assertTrue(NmfRun.signaturesEqual(sigs1, sigs2));
    }

}
