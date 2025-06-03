package com.hartwig.hmftools.purple.fitting;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.purple.fittingsnv.PeakModelData;
import com.hartwig.hmftools.purple.somatic.SubclonalLikelihood;

import org.junit.Test;

public class SubclonalLikelihoodFactoryTest
{
    @Test
    public void testSubclonalLikelihood()
    {
        List<PeakModelData> model = Lists.newArrayList(
                createModel(true, true, 0.0, 1),
                createModel(true, true, 0.05, 2),
                createModel(false, true, 0.05, 1),
                createModel(false, true, 0.05, 1),
                createModel(false, true, 0.1, 1));

        SubclonalLikelihood subclonalLikelihood = new SubclonalLikelihood(0.05, model);
        assertEquals(1, subclonalLikelihood.subclonalLikelihood(0.01), 0.01);

        assertEquals(0.5, subclonalLikelihood.subclonalLikelihood(0.025), 0.01);
        assertEquals(0.5, subclonalLikelihood.subclonalLikelihood(0.05), 0.01);
        assertEquals(0.5, subclonalLikelihood.subclonalLikelihood(0.0749), 0.01);
    }

    @Test
    public void testIgnoreInvalid()
    {
        List<PeakModelData> model = Lists.newArrayList(
                createModel(true, true, 0.05, 2),
                createModel(false, false, 0.05, 1),
                createModel(false, true, 0.05, 1));

        SubclonalLikelihood subclonalLikelihood = new SubclonalLikelihood(0.05, model);
        assertEquals(2d / 3d, subclonalLikelihood.subclonalLikelihood(0.05), 0.01);
    }

    @Test
    public void testOutsideRange()
    {
        SubclonalLikelihood subclonalLikelihood = new SubclonalLikelihood(0.05, Lists.newArrayList());
        assertEquals(0, subclonalLikelihood.subclonalLikelihood(0.05), 0.01);
    }

    private static PeakModelData createModel(boolean subclonal, boolean isValid, double bucket, double bucketWeight)
    {
        return new PeakModelData(0, 0, bucket, bucketWeight, isValid, subclonal);
    }
}
