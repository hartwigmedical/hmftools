package com.hartwig.hmftools.wisp.variant;

import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_CONFIG;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_REFERENCE_ID;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_RESULTS_WRITER;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_SAMPLE;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_SAMPLE_ID;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_TUMOR_ID;
import static com.hartwig.hmftools.wisp.variant.MosaicRoutines.findOutlierVariants;
import static com.hartwig.hmftools.wisp.variant.SomaticVariantUtils.addBqrData;
import static com.hartwig.hmftools.wisp.variant.SomaticVariantUtils.createSomaticVariant;
import static com.hartwig.hmftools.wisp.variant.SomaticVariantUtils.createVariantContext;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.wisp.WispConfig;

import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;

public class MosaicVariantTest
{
    @Test
    public void testOutlierVariants()
    {
        String ref = "A";
        String alt = "C";

        // add the non-outlier variants
        List<SomaticVariant> filteredVariants = Lists.newArrayList();

        int tumorDp = 100;
        int tumorAd = 50;

        for(int i = 0; i < 5; ++i)
        {
            filteredVariants.add(createSomaticVariant(ref, alt, tumorDp, tumorAd, 1000, 1));
        }

        for(int i = 0; i < 100; ++i)
        {
            filteredVariants.add(createSomaticVariant(ref, alt, tumorDp, tumorAd, 1000, 0));
        }

        SomaticVariant var1 = createSomaticVariant(ref, alt, tumorDp, tumorAd, 1000, 450);

        filteredVariants.add(var1);

        int totalSampleAd = filteredVariants.stream().mapToInt(x -> x.findSampleFragmentData(TEST_SAMPLE_ID).AlleleCount).sum();

        List<SomaticVariant> outliers = findOutlierVariants(TEST_SAMPLE_ID, filteredVariants, totalSampleAd);
        assertEquals(1, outliers.size());
        assertEquals(105, filteredVariants.size());
        assertTrue(outliers.contains(var1));

        // test again with a second outlier requiring use of the lower threshold
        filteredVariants.add(var1);

        SomaticVariant var2 = createSomaticVariant(ref, alt, tumorDp, tumorAd, 1000, 180);

        filteredVariants.add(var2);

        totalSampleAd = filteredVariants.stream().mapToInt(x -> x.findSampleFragmentData(TEST_SAMPLE_ID).AlleleCount).sum();

        outliers = findOutlierVariants(TEST_SAMPLE_ID, filteredVariants, totalSampleAd);
        assertEquals(2, outliers.size());
        assertEquals(105, filteredVariants.size());

        assertTrue(outliers.contains(var1));
        assertTrue(outliers.contains(var2));
        assertFalse(filteredVariants.contains(var1));
        assertFalse(filteredVariants.contains(var2));
    }

    @Test
    public void testMosaicVariants()
    {
        SomaticVariants somaticVariants = new SomaticVariants(TEST_CONFIG, TEST_RESULTS_WRITER, TEST_SAMPLE);

        String alt = "A";
        String ref = "C";
        addBqrData(somaticVariants.bqrAdjustment(), alt, 100_000, 10, 0);

        List<String> sampleIds = List.of(TEST_TUMOR_ID, TEST_REFERENCE_ID, TEST_SAMPLE_ID);

        // build specific variants
        int tumorDp = 100;
        int tumorAd = 50;
        double copyNumber = 2;
        int sampleDualDp = 400;
        int germlineAd = 2; // to fail germline filter

        VariantContext varContext1 = createVariantContext(
                ref, alt, 30, germlineAd, tumorDp, tumorAd, 1000, 0, sampleDualDp,
                0, copyNumber, copyNumber);

        somaticVariants.processVariant(sampleIds, varContext1);

        VariantContext varContext2 = createVariantContext(
                ref, alt, 30, 0, tumorDp, tumorAd, 1000, 1, sampleDualDp,
                0, copyNumber, copyNumber);

        varContext2.getCommonInfo().putAttribute(VariantTier.TIER, VariantTier.LOW_CONFIDENCE); // a filter that still counts for mosaics
        somaticVariants.processVariant(sampleIds, varContext2);

        VariantContext varContext3 = createVariantContext(
                ref, alt, 30, germlineAd, tumorDp, tumorAd, 1000, 0, sampleDualDp,
                0, copyNumber, copyNumber);

        somaticVariants.processVariant(sampleIds, varContext3);

        VariantContext varContext4 = createVariantContext(
                ref, alt, 30, 0, tumorDp, tumorAd, 1000, 1, sampleDualDp,
                0, copyNumber, copyNumber);

        varContext4.getCommonInfo().putAttribute(VariantTier.TIER, VariantTier.LOW_CONFIDENCE);
        somaticVariants.processVariant(sampleIds, varContext4);

        VariantContext varContext5 = createVariantContext(
                ref, alt, 30, 0, tumorDp, tumorAd, 1000, 0, sampleDualDp,
                0, copyNumber, copyNumber);

        varContext5.getCommonInfo().putAttribute(VariantTier.TIER, VariantTier.LOW_CONFIDENCE);
        somaticVariants.processVariant(sampleIds, varContext5);

        copyNumber = 2.5;
        VariantContext varContext6 = createVariantContext(
                ref, alt, 30, 0, tumorDp, tumorAd, 1000, 25, sampleDualDp,
                3, copyNumber, copyNumber);

        somaticVariants.processVariant(sampleIds, varContext6);

        VariantContext varContext7 = createVariantContext(
                ref, alt, 30, germlineAd, tumorDp, tumorAd, 1500, 2, sampleDualDp,
                1, copyNumber, copyNumber);

        varContext7.getCommonInfo().putAttribute(VariantTier.TIER, VariantTier.LOW_CONFIDENCE);
        somaticVariants.processVariant(sampleIds, varContext7);

        VariantContext varContext8 = createVariantContext(
                ref, alt, 30, 0, tumorDp, tumorAd, 1500, 0, sampleDualDp,
                0, copyNumber, copyNumber);

        somaticVariants.processVariant(sampleIds, varContext8);

        VariantContext varContext9 = createVariantContext(
                ref, alt, 30, 0, tumorDp, tumorAd, 1500, 1, sampleDualDp,
                1, copyNumber, copyNumber);

        somaticVariants.processVariant(sampleIds, varContext9);

        VariantContext varContext10 = createVariantContext(
                ref, alt, 30, 0, tumorDp, tumorAd, 1500, 0, sampleDualDp,
                0, copyNumber, copyNumber);

        somaticVariants.processVariant(sampleIds, varContext10);

        copyNumber = 3;
        VariantContext varContext11 = createVariantContext(
                ref, alt, 30, 0, tumorDp, tumorAd, 2000, 0, sampleDualDp,
                0, copyNumber, copyNumber);

        somaticVariants.processVariant(sampleIds, varContext11);

        VariantContext varContext12 = createVariantContext(
                ref, alt, 30, 0, tumorDp, tumorAd, 2000, 1, sampleDualDp,
                1, copyNumber, copyNumber);

        somaticVariants.processVariant(sampleIds, varContext12);

        VariantContext varContext13 = createVariantContext(
                ref, alt, 30, 0, tumorDp, tumorAd, 2000, 2, sampleDualDp,
                1, 15, 15);

        somaticVariants.processVariant(sampleIds, varContext13);

        VariantContext varContext14 = createVariantContext(
                ref, alt, 30, 0, tumorDp, tumorAd, 2000, 0, sampleDualDp,
                0, copyNumber, copyNumber);

        somaticVariants.processVariant(sampleIds, varContext14);

        VariantContext varContext15 = createVariantContext(
                ref, alt, 30, 0, tumorDp, tumorAd, 2000, 1, sampleDualDp,
                0, copyNumber, copyNumber);

        WispConfig.RunMosaicLogic = true;

        somaticVariants.processVariant(sampleIds, varContext15);

        SomaticPurityResult purityResult = somaticVariants.processSample(TEST_SAMPLE_ID);

        assertEquals(0.00396, purityResult.PurityCalcs.ErrorRate, 0.0001);
        assertEquals(56, purityResult.mosaicFragments(), 1);
        assertEquals(12, purityResult.mosaicDualFragments(), 1);

        WispConfig.RunMosaicLogic = false;

    }
}
