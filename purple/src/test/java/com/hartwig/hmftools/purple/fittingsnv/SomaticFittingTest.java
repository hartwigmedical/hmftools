package com.hartwig.hmftools.purple.fittingsnv;

import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_FIT_TUMOR_ONLY_HOTSPOT_VAF_CUTOFF;
import static com.hartwig.hmftools.purple.fittingsnv.SomaticPurityFitter.calc75thPercentileValue;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;

import org.junit.Test;

public class SomaticFittingTest
{
    @Test
    public void testSomaticVafSelection()
    {
        List<Double> vafs = Lists.newArrayList(0.0, 0.25, 0.5, 0.75, 1.0);

        double calcVaf = calc75thPercentileValue(vafs);
        assertEquals(0.75, calcVaf, 0.01);

        vafs = Lists.newArrayList(0.0, 0.33, 0.66, 1.0);

        calcVaf = calc75thPercentileValue(vafs);
        assertEquals(0.75, calcVaf, 0.01);

        vafs.clear();
        for(double vaf = 0; vaf <= 0.5; vaf = vaf + 0.05)
        {
            vafs.add(vaf);
        }

        calcVaf = calc75thPercentileValue(vafs);
        assertEquals(0.375, calcVaf, 0.01);

        vafs = Lists.newArrayList(0.25, SOMATIC_FIT_TUMOR_ONLY_HOTSPOT_VAF_CUTOFF);
        calcVaf = calc75thPercentileValue(vafs);
        assertEquals(0.25, calcVaf, 0.01);

        double topVaf = SOMATIC_FIT_TUMOR_ONLY_HOTSPOT_VAF_CUTOFF - 0.01;
        vafs = Lists.newArrayList(0.25, topVaf);
        calcVaf = calc75thPercentileValue(vafs);
        assertEquals(topVaf, calcVaf, 0.01);

        vafs = Lists.newArrayList(0.5);
        calcVaf = calc75thPercentileValue(vafs);
        assertEquals(0.5, calcVaf, 0.01);
    }
}
