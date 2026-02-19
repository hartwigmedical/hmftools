package com.hartwig.hmftools.amber.contamination;

import com.hartwig.hmftools.common.amber.BaseDepthData;
import com.hartwig.hmftools.common.amber.ImmutableBaseDepthData;

import org.junit.Assert;
import org.junit.Test;

public class TumorContaminationTest
{
    @Test
    public void tumorVafTest()
    {
        BaseDepthData tumorBdd = ImmutableBaseDepthData.builder()
                .refSupport(100)
                .alt(BaseDepthData.Base.T)
                .ref(BaseDepthData.Base.G)
                .altSupport(5)
                .readDepth(100)
                .indelCount(12)
                .build();
        TumorContamination tc = new TumorContamination("chr1", 1000, null, tumorBdd);
        Assert.assertEquals(0.05, tc.tumorVaf(), 0.001);
    }
}
