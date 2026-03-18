package com.hartwig.hmftools.amber.contamination;

import com.hartwig.hmftools.common.amber.AmberBase;
import com.hartwig.hmftools.common.amber.BaseDepthData;

import org.junit.Assert;
import org.junit.Test;

public class TumorContaminationTest
{
    @Test
    public void tumorVafTest()
    {
        BaseDepthData tumorBdd = new BaseDepthData(
                AmberBase.T,
                AmberBase.G,
                100,
                0,
                95,
                5);
        TumorContamination tc = new TumorContamination("chr1", 1000, null, tumorBdd);
        Assert.assertEquals(0.05, tc.tumorVaf(), 0.001);
    }
}
