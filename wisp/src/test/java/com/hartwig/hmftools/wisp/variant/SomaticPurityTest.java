package com.hartwig.hmftools.wisp.variant;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_VARIANT_CN;
import static com.hartwig.hmftools.common.variant.SageVcfTags.AVG_EDGE_DISTANCE_PERC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.AVG_RECALIBRATED_BASE_QUAL;
import static com.hartwig.hmftools.common.variant.SageVcfTags.NEARBY_INDEL_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_STRAND_BIAS;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TRINUCLEOTIDE_CONTEXT;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_CONFIG;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_REFERENCE_ID;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_RESULTS_WRITER;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_SAMPLE;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_SAMPLE_ID;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_TUMOR_ID;
import static com.hartwig.hmftools.wisp.variant.ClonalityMethod.NONE;
import static com.hartwig.hmftools.wisp.variant.FilterReason.AVG_EDGE_DIST;
import static com.hartwig.hmftools.wisp.variant.FilterReason.DUAL_ERROR_RATE;
import static com.hartwig.hmftools.wisp.variant.FilterReason.GC_RATIO;
import static com.hartwig.hmftools.wisp.variant.FilterReason.GERMLINE_AF;
import static com.hartwig.hmftools.wisp.variant.FilterReason.LOW_CONFIDENCE;
import static com.hartwig.hmftools.wisp.variant.FilterReason.LOW_QUAL_PER_AD;
import static com.hartwig.hmftools.wisp.variant.FilterReason.MAPPABILITY;
import static com.hartwig.hmftools.wisp.variant.FilterReason.NEARBY_INDEL;
import static com.hartwig.hmftools.wisp.variant.FilterReason.NON_SNV;
import static com.hartwig.hmftools.wisp.variant.FilterReason.NO_PASS;
import static com.hartwig.hmftools.wisp.variant.FilterReason.REPEAT_COUNT;
import static com.hartwig.hmftools.wisp.variant.FilterReason.SUBCLONAL;
import static com.hartwig.hmftools.wisp.variant.SomaticVariantUtils.DEFAULT_TUMOR_DEPTH;
import static com.hartwig.hmftools.wisp.variant.SomaticVariantUtils.addBqrData;
import static com.hartwig.hmftools.wisp.variant.SomaticVariantUtils.createSomaticVariant;
import static com.hartwig.hmftools.wisp.variant.SomaticVariantUtils.createVariantContext;
import static com.hartwig.hmftools.wisp.variant.SomaticVariantUtils.formTriNucContext;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.variant.PaveVcfTags;
import com.hartwig.hmftools.common.variant.SageVcfTags;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.wisp.WispConfig;

import org.junit.Test;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;

public class SomaticPurityTest
{
    private final BqrAdjustment mBbqrAdjustment;
    private final SomaticPurityEstimator mPurityEstimator;

    public SomaticPurityTest()
    {
        mBbqrAdjustment = new BqrAdjustment(TEST_CONFIG);

        mPurityEstimator = new SomaticPurityEstimator(TEST_CONFIG, TEST_RESULTS_WRITER, TEST_SAMPLE, mBbqrAdjustment);

        addBqrData(mBbqrAdjustment, "C", 100_000, 10, 0);
    }

    @Test
    public void testVariantFilters()
    {
        Map<String,Object> commonAttributes = Maps.newHashMap();
        Map<String,Object> refAttributes = Maps.newHashMap();
        Map<String,Object> tumorAttributes = Maps.newHashMap();
        Map<String,Object> sampleAttributes = Maps.newHashMap();

        refAttributes.put(AVG_RECALIBRATED_BASE_QUAL, "0,0" );

        String ref = "A";
        String alt = "C";

        // first with filter check passing
        commonAttributes.put(TRINUCLEOTIDE_CONTEXT, formTriNucContext(alt));
        commonAttributes.put(PURPLE_CN, 2.0);
        commonAttributes.put(PURPLE_VARIANT_CN, 1.0);

        sampleAttributes.put(AVG_EDGE_DISTANCE_PERC, "75,75");
        sampleAttributes.put(READ_STRAND_BIAS, "0.5,0.5");

        WispConfig config = new WispConfig(new MockRefGenome(), 0.4);
        SomaticVariants somaticVariants = new SomaticVariants(config, TEST_RESULTS_WRITER, TEST_SAMPLE);

        VariantContext variantContext = createVariantContext(
                ref, alt, 30, 0, 100, 50, 1000, 100, 100, 5,
                commonAttributes, refAttributes, tumorAttributes, sampleAttributes);

        VariantContextDecorator variant = new VariantContextDecorator(variantContext);

        List<FilterReason> filters = somaticVariants.checkFilters(variant, 0, 0.5);

        assertTrue(filters.isEmpty());

        // now trigger each filter
        commonAttributes.put(VariantTier.TIER, VariantTier.LOW_CONFIDENCE);
        commonAttributes.put(PaveVcfTags.MAPPABILITY, 0.8);
        commonAttributes.put(NEARBY_INDEL_FLAG, true);
        commonAttributes.put(SageVcfTags.REPEAT_COUNT, 10);
        commonAttributes.put(PURPLE_VARIANT_CN, 0.5);

        alt = "AC";
        variantContext = createVariantContext(
                ref, alt, 30, 2, 100, 50, 1000, 100, 100, 5,
                commonAttributes, refAttributes, tumorAttributes, sampleAttributes);

        variantContext.getCommonInfo().addFilter("SomeFilter");

        variant = new VariantContextDecorator(variantContext);
        filters = somaticVariants.checkFilters(variant, 0.7, 0.1);

        assertTrue(filters.contains(GERMLINE_AF));
        assertTrue(filters.contains(NEARBY_INDEL));
        assertTrue(filters.contains(NO_PASS));
        assertTrue(filters.contains(REPEAT_COUNT));
        assertTrue(filters.contains(NON_SNV));
        assertTrue(filters.contains(MAPPABILITY));
        assertTrue(filters.contains(LOW_CONFIDENCE));
        assertTrue(filters.contains(SUBCLONAL));
        assertTrue(filters.contains(GC_RATIO));
    }

    @Test
    public void testSampleFilters()
    {
        Map<String,Object> commonAttributes = Maps.newHashMap();
        Map<String,Object> refAttributes = Maps.newHashMap();
        Map<String,Object> tumorAttributes = Maps.newHashMap();
        Map<String,Object> sampleAttributes = Maps.newHashMap();

        String ref = "A";
        String alt = "C";

        // first with filter check passing
        commonAttributes.put(TRINUCLEOTIDE_CONTEXT, formTriNucContext(alt));
        commonAttributes.put(PURPLE_CN, 2.0);
        commonAttributes.put(PURPLE_VARIANT_CN, 1.0);

        sampleAttributes.put(AVG_EDGE_DISTANCE_PERC, "75,75");
        sampleAttributes.put(READ_STRAND_BIAS, "0.5,0.5");

        WispConfig config = new WispConfig(new MockRefGenome(), 0.4);
        SomaticVariants somaticVariants = new SomaticVariants(config, TEST_RESULTS_WRITER, TEST_SAMPLE);

        SomaticVariant variant = createSomaticVariant(ref, alt,100, 50, 100, 0);

        int sampleDp = 100;
        int sampleAd = 10;
        UmiTypeCounts umiTypeCounts = new UmiTypeCounts(100, 0, 0, 10, 0, 0);

        String sampleId = "SAMPLE_002";

        Genotype gtSample = new GenotypeBuilder().attributes(sampleAttributes).name(sampleId).DP(-1)
                .AD(new int[] { sampleDp, sampleDp - sampleAd, sampleAd }).noPL().GQ(-1).make();

        GenotypeFragments sampleFragData = new GenotypeFragments(
                sampleId, 10, 100, 370, umiTypeCounts, gtSample);

        somaticVariants.checkSampleDataFilters(variant, sampleFragData);

        assertTrue(sampleFragData.filterReasons().isEmpty());

        sampleAttributes.put(AVG_EDGE_DISTANCE_PERC, "0.05,0.05");
        // sampleAttributes.put(READ_STRAND_BIAS, "0.5,0.5");

        gtSample = new GenotypeBuilder().attributes(sampleAttributes).name(sampleId).DP(-1)
                .AD(new int[] { sampleDp, sampleDp - sampleAd, sampleAd }).noPL().GQ(-1).make();

        sampleFragData = new GenotypeFragments(
                sampleId, 10, 100, 200, umiTypeCounts, gtSample);

        addBqrData(somaticVariants.bqrAdjustment(), alt, 10_000, 10, 10);

        somaticVariants.checkSampleDataFilters(variant, sampleFragData);

        assertTrue(sampleFragData.filterReasons().contains(LOW_QUAL_PER_AD));
        assertTrue(sampleFragData.filterReasons().contains(AVG_EDGE_DIST));
        assertTrue(sampleFragData.dualFilterReasons().contains(DUAL_ERROR_RATE));
    }

    @Test
    public void testSomaticVariantsBasic()
    {
        SomaticVariants somaticVariants = new SomaticVariants(TEST_CONFIG, TEST_RESULTS_WRITER, TEST_SAMPLE);

        String alt = "C";
        addBqrData(somaticVariants.bqrAdjustment(), alt, 100_000, 10, 0);

        List<String> sampleIds = List.of(TEST_TUMOR_ID, TEST_REFERENCE_ID, TEST_SAMPLE_ID);

        for(int i = 0; i < 10; ++i)
        {
            VariantContext variantContext = createVariantContext(
                    "A", "C", DEFAULT_TUMOR_DEPTH, 50, 200, 100);

            somaticVariants.processVariant(sampleIds, variantContext);
        }

        SomaticPurityResult purityResult = somaticVariants.processSample(TEST_SAMPLE_ID);

        assertEquals(10, purityResult.TotalVariants);
        assertEquals(SnvFitStatus.OK, purityResult.Status);
        assertTrue(purityResult.CapMaxTumorCopyNumber);

        assertEquals(0.50, purityResult.PurityCalcs.RawPurityEstimate, 0.01);
        assertEquals(0.54, purityResult.PurityCalcs.PurityRangeHigh, 0.01);
        assertEquals(0.46, purityResult.PurityCalcs.PurityRangeLow, 0.01);
        assertEquals(0, purityResult.PurityCalcs.Probability, 0.01);
        assertEquals(0.00097, purityResult.PurityCalcs.LodPurityEstimate, 0.001);
        assertEquals(0, purityResult.PurityCalcs.BqrQualThreshold);
        assertEquals(NONE, purityResult.PurityCalcs.Clonality.Method);
        assertEquals(3.2e-5, purityResult.PurityCalcs.ErrorRate, 1e-4);
    }

    @Test
    public void testSomaticVariantsFilters()
    {
        SomaticVariants somaticVariants = new SomaticVariants(TEST_CONFIG, TEST_RESULTS_WRITER, TEST_SAMPLE);

    }

   @Test
    public void testStandardSomaticFit()
    {
        // mPurityEstimator.calculatePurity(TEST_SAMPLE_ID, )

    }
}
