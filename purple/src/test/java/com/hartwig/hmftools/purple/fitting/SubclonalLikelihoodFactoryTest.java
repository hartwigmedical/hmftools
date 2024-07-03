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
        final List<PeakModelData> model = Lists.newArrayList(
                createModel(true, true, 0.0, 1),
                createModel(true, true, 0.05, 2),
                createModel(false, true, 0.05, 1),
                createModel(false, true, 0.05, 1),
                createModel(false, true, 0.1, 1));

        final SubclonalLikelihood victim = new SubclonalLikelihood(0.05, model);
        assertEquals(1, victim.subclonalLikelihood(0.01), 0.01);

        assertEquals(0.5, victim.subclonalLikelihood(0.025), 0.01);
        assertEquals(0.5, victim.subclonalLikelihood(0.05), 0.01);
        assertEquals(0.5, victim.subclonalLikelihood(0.0749), 0.01);
    }

    @Test
    public void testIgnoreInvalid()
    {
        final List<PeakModelData> model = Lists.newArrayList(
                createModel(true, true, 0.05, 2),
                createModel(false, false, 0.05, 1),
                createModel(false, true, 0.05, 1));

        final SubclonalLikelihood victim = new SubclonalLikelihood(0.05, model);
        assertEquals(2d / 3d, victim.subclonalLikelihood(0.05), 0.01);
    }

    @Test
    public void testOutsideRange()
    {
        final SubclonalLikelihood victim = new SubclonalLikelihood(0.05, Lists.newArrayList());
        assertEquals(0, victim.subclonalLikelihood(0.05), 0.01);
    }

    private static PeakModelData createModel(boolean subclonal, boolean isValid, double bucket, double bucketWeight)
    {
        return new PeakModelData(0, 0, bucket, bucketWeight, isValid, subclonal);
    }
}
