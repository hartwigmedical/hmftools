package com.hartwig.hmftools.wisp.variant;

import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_CONFIG;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_RESULTS_WRITER;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_SAMPLE;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_SAMPLE_ID;
import static com.hartwig.hmftools.wisp.variant.SomaticVariantUtils.addBqrData;
import static com.hartwig.hmftools.wisp.variant.SomaticVariantUtils.createSomaticVariant;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;

import org.junit.Test;

public class ClonalModelsTest
{
    @Test
    public void vafPeakModelTest()
    {
        String ref = "A";
        String alt = "C";

        // add the non-outlier variants
        List<SomaticVariant> variants = Lists.newArrayList();

        int tumorDp = 100;
        int tumorAd = 50;
        int sampleDp = 1000;

        for(int fragCount = 1; fragCount <= 5; ++fragCount)
        {
            for(int i = 0; i < 10; ++i)
            {
                variants.add(createSomaticVariant(ref, alt, tumorDp, tumorAd, sampleDp, fragCount));
            }
        }

        MosaicData mosaicData = new MosaicData(0, 0, 0);

        BqrAdjustment bqrAdjustment = new BqrAdjustment(TEST_CONFIG);

        addBqrData(bqrAdjustment, alt, 100_000, 10, 0);

        SomaticPurityEstimator estimator = new SomaticPurityEstimator(TEST_CONFIG, TEST_RESULTS_WRITER, TEST_SAMPLE, bqrAdjustment);

        SomaticPurityResult purityResult = estimator.calculatePurity(TEST_SAMPLE_ID, variants, variants.size(), mosaicData, false);

        assertEquals(ClonalityMethod.VAF_PEAK, purityResult.PurityCalcs.Clonality.Method);

        assertEquals(0.003, purityResult.PurityCalcs.Clonality.Vaf, 0.001);
        assertEquals(50, purityResult.PurityCalcs.Clonality.VariantCount);
        assertEquals(0.003, purityResult.PurityCalcs.PurityEstimate, 0.001);
        assertEquals(SnvFitStatus.OK, purityResult.Status);
    }

    @Test
    public void lowCountModelTest()
    {
        String ref = "A";
        String alt = "C";

        // add the non-outlier variants
        List<SomaticVariant> variants = Lists.newArrayList();

        int tumorDp = 100;
        int tumorAd = 50;

        int sampleDp = 20;

        for(int i = 0; i < 10; ++i)
        {
            variants.add(createSomaticVariant(ref, alt, tumorDp, tumorAd, sampleDp, 0));
        }

        for(int i = 0; i < 4; ++i)
        {
            variants.add(createSomaticVariant(ref, alt, tumorDp, tumorAd, sampleDp, 1));
        }

        for(int i = 0; i < 5; ++i)
        {
            variants.add(createSomaticVariant(ref, alt, tumorDp, tumorAd, sampleDp, 2));
        }

        for(int i = 0; i < 3; ++i)
        {
            variants.add(createSomaticVariant(ref, alt, tumorDp, tumorAd, sampleDp, 3));
        }

        MosaicData mosaicData = new MosaicData(0, 0, 0);

        BqrAdjustment bqrAdjustment = new BqrAdjustment(TEST_CONFIG);

        addBqrData(bqrAdjustment, alt, 100_000, 10, 0);

        SomaticPurityEstimator estimator = new SomaticPurityEstimator(TEST_CONFIG, TEST_RESULTS_WRITER, TEST_SAMPLE, bqrAdjustment);

        SomaticPurityResult purityResult = estimator.calculatePurity(TEST_SAMPLE_ID, variants, variants.size(), mosaicData, false);

        assertEquals(ClonalityMethod.LOW_COUNT, purityResult.PurityCalcs.Clonality.Method);

        assertEquals(0.0936, purityResult.PurityCalcs.Clonality.Vaf, 0.001);
        assertEquals(12, purityResult.PurityCalcs.Clonality.VariantCount);
        assertEquals(0.0936, purityResult.PurityCalcs.PurityEstimate, 0.001);
        assertEquals(SnvFitStatus.OK, purityResult.Status);
    }
}
