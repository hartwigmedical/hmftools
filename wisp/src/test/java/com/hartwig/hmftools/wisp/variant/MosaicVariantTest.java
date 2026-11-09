package com.hartwig.hmftools.wisp.variant;

import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_CONFIG;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_REFERENCE_ID;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_RESULTS_WRITER;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_SAMPLE;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_SAMPLE_ID;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_TUMOR_ID;
import static com.hartwig.hmftools.wisp.variant.SomaticVariantUtils.addBqrData;
import static com.hartwig.hmftools.wisp.variant.SomaticVariantUtils.createVariantContext;

import java.util.List;

import com.hartwig.hmftools.common.variant.VariantTier;

import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;

public class MosaicVariantTest
{
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
                0, copyNumber, copyNumber);

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
                0, copyNumber, copyNumber);

        somaticVariants.processVariant(sampleIds, varContext12);

        VariantContext varContext13 = createVariantContext(
                ref, alt, 30, 0, tumorDp, tumorAd, 2000, 2, sampleDualDp,
                0, 15, 15);

        somaticVariants.processVariant(sampleIds, varContext13);

        VariantContext varContext14 = createVariantContext(
                ref, alt, 30, 0, tumorDp, tumorAd, 2000, 0, sampleDualDp,
                0, copyNumber, copyNumber);

        somaticVariants.processVariant(sampleIds, varContext14);

        VariantContext varContext15 = createVariantContext(
                ref, alt, 30, 0, tumorDp, tumorAd, 2000, 1, sampleDualDp,
                0, copyNumber, copyNumber);

        somaticVariants.processVariant(sampleIds, varContext15);

        SomaticPurityResult purityResult = somaticVariants.processSample(TEST_SAMPLE_ID);

        //  assertEquals(10, purityResult.TotalVariants);



    }
}
