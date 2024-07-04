package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_VARIANT_CN;
import static com.hartwig.hmftools.purple.MiscTestUtils.REF_SAMPLE_ID;
import static com.hartwig.hmftools.purple.MiscTestUtils.SAMPLE_ID;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.SegmentSupport;

import org.junit.Test;

import java.util.List;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import junit.framework.TestCase;

public class SomaticPurityEnrichmentTest extends TestCase
{
    @Test
    public void testCalculateBiallelic()
    {
        // 27.06.2024: New unit tests due to model changes
        
        // Scenario 1:
        // input values:
        double CN = 5.1;
        double MACN = 2;
        double VCN = 2.83;
        int alleleReadCount = 42;

        // expected output:
        double expectedBiallelicProbability = 0.02;

        // MACN = (1 - BAF) * CN <=> BAF = 1 - MACN / CN
        double correspondingBAF = 1 - (MACN / CN);

        PurpleCopyNumber copyNumber = createCopyNumber(CN, correspondingBAF);

        SomaticVariant variant = createVariant(alleleReadCount, VCN);

        double biallelicProbability = SomaticPurityEnrichment.calculateBiallelic(copyNumber, variant);
        assertEquals(expectedBiallelicProbability, biallelicProbability, 0.01);

        // Scenario 2:
        // input values:
        CN = 2.01;
        MACN = 0.0004;
        VCN = 2.11;
        alleleReadCount = 52;

        // expected output:
        expectedBiallelicProbability = 0.995;
        
        correspondingBAF = 1 - (MACN / CN);
        copyNumber = createCopyNumber(CN, correspondingBAF);
        variant = createVariant(alleleReadCount, VCN);

        biallelicProbability = SomaticPurityEnrichment.calculateBiallelic(copyNumber, variant);
        assertEquals(expectedBiallelicProbability, biallelicProbability, 0.01);

        // Scenario 3
        // input values:
        CN = 2.98;
        MACN = 0.991;
        VCN = 2.33;
        alleleReadCount = 36;

        // expected output:
        expectedBiallelicProbability = 0.054;

        correspondingBAF = 1 - (MACN / CN);
        copyNumber = createCopyNumber(CN, correspondingBAF);
        variant = createVariant(alleleReadCount, VCN);

        biallelicProbability = SomaticPurityEnrichment.calculateBiallelic(copyNumber, variant);
        assertEquals(expectedBiallelicProbability, biallelicProbability, 0.01);

        // Scenario 4
        // input values:
        CN = 43.78;
        MACN = 0.0;
        VCN = 42.74;
        alleleReadCount = 627;

        // expected output:
        expectedBiallelicProbability = 0.442;

        correspondingBAF = 1 - (MACN / CN);
        copyNumber = createCopyNumber(CN, correspondingBAF);
        variant = createVariant(alleleReadCount, VCN);

        biallelicProbability = SomaticPurityEnrichment.calculateBiallelic(copyNumber, variant);
        assertEquals(expectedBiallelicProbability, biallelicProbability, 0.01);

        // Scenario 5 - initially called BIALLELIC in COLO829v003T chr1:50084056
        // input values:
        CN = 2.81;
        MACN = 0.811;
        VCN = 2.83;
        alleleReadCount = 73;

        // expected output:
        expectedBiallelicProbability = 0.833;

        correspondingBAF = 1 - (MACN / CN);
        copyNumber = createCopyNumber(CN, correspondingBAF);
        variant = createVariant(alleleReadCount, VCN);

        biallelicProbability = SomaticPurityEnrichment.calculateBiallelic(copyNumber, variant);
        assertEquals(expectedBiallelicProbability, biallelicProbability, 0.01);
    }
    
    private static SomaticVariant createVariant(int alleleReadCount, double variantCopyNumber)
    {
        GenotypeBuilder genotypeBuilder = new GenotypeBuilder(SAMPLE_ID);
        genotypeBuilder.AD(new int[] { alleleReadCount, alleleReadCount });
        genotypeBuilder.alleles(Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL));
        Genotype genotype = genotypeBuilder.make();
        
        List<Allele> alleles = List.of(Allele.create("A", true), Allele.create("T", false));
        
        VariantContextBuilder builder = new VariantContextBuilder();
        builder.chr(CHR_1).start(100);
        builder.alleles(alleles);
        builder.computeEndFromAlleles(alleles, 100);

        builder.genotypes(List.of(genotype));

        builder.attribute(PURPLE_VARIANT_CN, variantCopyNumber);

        VariantContext variantContext = builder.make();

        return new SomaticVariant(variantContext, SAMPLE_ID, REF_SAMPLE_ID);
    }

    private static PurpleCopyNumber createCopyNumber(final double copyNumber, final double averageActualBAF)
    {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome("chr1")
                .start(1)
                .end(1000)
                .averageTumorCopyNumber(copyNumber)
                .segmentStartSupport(SegmentSupport.NONE)
                .segmentEndSupport(SegmentSupport.NONE)
                .method(CopyNumberMethod.UNKNOWN)
                .bafCount(0)
                .depthWindowCount(1)
                .gcContent(0)
                .minStart(1)
                .maxStart(10)
                .averageObservedBAF(0.5)
                .averageActualBAF(averageActualBAF)
                .build();
    }
}