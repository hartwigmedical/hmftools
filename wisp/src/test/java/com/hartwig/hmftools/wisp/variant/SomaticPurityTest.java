package com.hartwig.hmftools.wisp.variant;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASES;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_VARIANT_CN;
import static com.hartwig.hmftools.common.variant.SageVcfTags.AVG_EDGE_DISTANCE_PERC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.AVG_RECALIBRATED_BASE_QUAL;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_QUALITY;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_STRAND_BIAS;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TRINUCLEOTIDE_CONTEXT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNTS;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_CONFIG;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_REFERENCE_ID;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_RESULTS_WRITER;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_SAMPLE;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_SAMPLE_ID;
import static com.hartwig.hmftools.wisp.purity.WispTestUtils.TEST_TUMOR_ID;
import static com.hartwig.hmftools.wisp.variant.ClonalityMethod.NONE;
import static com.hartwig.hmftools.wisp.variant.GenotypeFragments.parseReadArrayAttribute;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.wisp.ResultsWriter;
import com.hartwig.hmftools.wisp.SampleData;
import com.hartwig.hmftools.wisp.WispConfig;
import com.hartwig.hmftools.wisp.variant.BqrAdjustment;
import com.hartwig.hmftools.wisp.variant.FilterReason;
import com.hartwig.hmftools.wisp.variant.SomaticPurityEstimator;
import com.hartwig.hmftools.wisp.variant.SomaticVariant;

import org.junit.Test;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

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
    public void tesSomaticVariantsBasic()
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

        assertEquals(1.00, purityResult.PurityCalcs.RawPurityEstimate, 0.01);
        assertEquals(1.00, purityResult.PurityCalcs.PurityRangeHigh, 0.01);
        assertEquals(0.92, purityResult.PurityCalcs.PurityRangeLow, 0.01);
        assertEquals(0, purityResult.PurityCalcs.Probability, 0.01);
        assertEquals(0.002, purityResult.PurityCalcs.LodPurityEstimate, 0.001);
        assertEquals(0, purityResult.PurityCalcs.BqrQualThreshold);
        assertEquals(NONE, purityResult.PurityCalcs.Clonality.Method);
        assertEquals(3.2e-5, purityResult.PurityCalcs.ErrorRate, 1e-4);
    }

    @Test
    public void testStandardSomaticFit()
    {
        // mPurityEstimator.calculatePurity(TEST_SAMPLE_ID, )

    }

    protected static final int DEFAULT_REF_DEPTH = 30;
    protected static final int DEFAULT_TUMOR_DEPTH = 100;
    protected static final int DEFAULT_BASE_QUAL = 37;

    private static void addReadCountAttributes(
            final Map<String,Object> attributes, int alleleFrags, int depth, int dualFrags, int dualDepth)
    {
        int tumorRefCount = depth - alleleFrags;

        attributes.put(READ_CONTEXT_COUNT, format("%d,0,0,0,%d,%d", alleleFrags, tumorRefCount, depth));

        attributes.put(READ_CONTEXT_QUALITY, format("%d,0,0,0,%d,%d",
                alleleFrags * DEFAULT_BASE_QUAL, tumorRefCount * DEFAULT_BASE_QUAL, depth * DEFAULT_BASE_QUAL));

        attributes.put(UMI_TYPE_COUNTS, format("%d,0,%d,%d,0,%d", depth - dualDepth, dualDepth, alleleFrags - dualFrags, dualFrags));
    }

    protected static String formTriNucContext(final String alt)
    {
        return "A" + alt + "T";
    }

    protected static VariantContext createVariantContext(
            final String ref, final String alt, int tumorDp, int tumorAd, int sampleDp, int sampleAd)
    {
        return createVariantContext(
                ref, alt, DEFAULT_REF_DEPTH, 0, tumorDp, tumorAd, sampleDp, sampleAd,
                0, 0, 2.0, 2.0);
    }

    protected static VariantContext createVariantContext(
            final String ref, final String alt, int refDp, int refAd, int tumorDp, int tumorAd,
            int sampleDp, int sampleAd, int sampleDualDp, int sampleDualAd,
            double copyNumber, double variantCopyNumber)
    {
        Map<String,Object> commonAttributes = Maps.newHashMap();
        Map<String,Object> refAttributes = Maps.newHashMap();
        Map<String,Object> tumorAttributes = Maps.newHashMap();
        Map<String,Object> sampleAttributes = Maps.newHashMap();

        refAttributes.put(AVG_RECALIBRATED_BASE_QUAL, "0,0" );

        commonAttributes.put(TRINUCLEOTIDE_CONTEXT, formTriNucContext(alt));
        commonAttributes.put(PURPLE_CN, copyNumber);
        commonAttributes.put(PURPLE_VARIANT_CN, variantCopyNumber);

        sampleAttributes.put(AVG_EDGE_DISTANCE_PERC, "75,75");
        sampleAttributes.put(READ_STRAND_BIAS, "0.5,0.5");

        return createVariantContext(
                ref, alt, refDp, refAd, tumorDp, tumorAd, sampleDp, sampleAd, sampleDualDp, sampleDualAd,
                commonAttributes, refAttributes, tumorAttributes, sampleAttributes);
    }

    protected static VariantContext createVariantContext(
            final String ref, final String alt, int refDp, int refAd, int tumorDp, int tumorAd,
            int sampleDp, int sampleAd, int sampleDualDp, int sampleDualAd,
            final Map<String,Object> commonAttributes, final Map<String,Object> refAttributes,
            final Map<String,Object> tumorAttributes, final Map<String,Object> sampleAttributes)
    {
        VariantContextBuilder builder = new VariantContextBuilder();

        List<Allele> alleles = Lists.newArrayList();
        alleles.add(Allele.create(ref, true));
        alleles.add(Allele.create(alt, false));

        refAd = min(refAd, refDp);
        int refRefCount = refDp - refAd;

        addReadCountAttributes(refAttributes, refAd, refDp, 0, 0);

        Genotype gtNormal = new GenotypeBuilder().attributes(refAttributes).name(TEST_REFERENCE_ID)
                .DP(refDp)
                .AD(new int[] { refRefCount, refAd })
                .noPL().GQ(-1).make();

        tumorAd = min(tumorAd, tumorDp);

        int tumorRefCount = tumorDp - tumorAd;
        addReadCountAttributes(tumorAttributes, tumorAd, tumorDp, 0, 0);

        Genotype gtTumor = new GenotypeBuilder().attributes(tumorAttributes).name(TEST_TUMOR_ID)
                .DP(DEFAULT_TUMOR_DEPTH)
                .AD(new int[] { tumorRefCount, tumorAd })
                .noPL().GQ(-1).make();

        sampleAd = min(sampleAd, sampleDp);
        int sampleRefCount = sampleDp - sampleAd;

        addReadCountAttributes(sampleAttributes, sampleAd, sampleDp, sampleDualAd, sampleDualDp);

        Genotype gtSample = new GenotypeBuilder().attributes(sampleAttributes).name(TEST_SAMPLE_ID)
                .DP(-1)
                .AD(new int[] { sampleDp, sampleRefCount, sampleAd })
                .noPL().GQ(-1).make();

        GenotypesContext genotypesContext = GenotypesContext.create(gtNormal, gtTumor, gtSample);

        double logError = -(1000 / 10.0);

        int start = 1000;
        int stop = start + ref.length() - 1;

        VariantContext variantContext = builder
                .source("SOURCE")
                .chr(CHR_1)
                .start(start)
                .stop(stop)
                .alleles(alleles)
                .genotypes(genotypesContext)
                .attributes(commonAttributes)
                .log10PError(logError)
                .unfiltered()
                .make(true);

        return variantContext;
    }

    protected SomaticVariant createSomaticVariant(
            final String ref, final String alt, int tumorDp, int tumorAd, int sampleDp, int sampleAd)
    {
        VariantContext variantContext = createVariantContext(ref, alt, tumorDp, tumorAd, sampleDp, sampleAd);

        VariantContextDecorator decorator = new VariantContextDecorator(variantContext);

        List<FilterReason> filters = Lists.newArrayList();

        return new SomaticVariant(decorator, 0, filters, false);
    }

    protected static void addBqrData(final BqrAdjustment bqrAdjustment, final String alt, int totalCount, int altCount, int altDualCount)
    {
        for(char base : DNA_BASES)
        {
            String baseStr = String.valueOf(base);
            String tnc = formTriNucContext(baseStr);
            BqrContextData bqrData = bqrAdjustment.getOrCreate(tnc, baseStr, ConsensusType.NONE);
            bqrData.TotalCount = totalCount;

            if(baseStr.equals(alt))
                bqrData.AltCount = altCount;

            if(altDualCount > 0)
            {
                bqrData = bqrAdjustment.getOrCreate(tnc, baseStr, ConsensusType.DUAL);
                bqrData.TotalCount = totalCount;

                if(baseStr.equals(alt))
                    bqrData.AltCount = altDualCount;
            }
        }
    }
}
