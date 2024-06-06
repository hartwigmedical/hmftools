package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_VARIANT_CN;
import static com.hartwig.hmftools.purple.TestUtils.REF_SAMPLE_ID;
import static com.hartwig.hmftools.purple.TestUtils.SAMPLE_ID;

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
        /* Scenario 1: High MACN -> No LOH -> expect low probability of biallelic */

        // input values:
        double CN = 3.08;
        double MACN = 1.01;
        double VCN = 0.941; // has totalReadCount as an input
        int alleleReadCount = 39;

        // expected output:
        double expectedBiallelicProbability = 0.01;

        // MACN = (1 - BAF) * CN <=> BAF = 1 - MACN / CN
        double correspondingBAF = 1 - (MACN / CN);

        PurpleCopyNumber copyNumber = createCopyNumber(CN, correspondingBAF);

        SomaticVariant variant = createVariant(alleleReadCount, VCN);

        double biallelicProbability = SomaticPurityEnrichment.calculateBiallelic(copyNumber, variant);
        assertEquals(expectedBiallelicProbability, biallelicProbability, 0.01);

        /* Scenario 2: Low MACN -> LOH -> expect high probability of biallelic */

        // input values:
        CN = 0.991;
        MACN = 0.0021;
        VCN = 1.06;
        alleleReadCount = 37;

        // expected output:
        expectedBiallelicProbability = 1.00;
        
        correspondingBAF = 1 - (MACN / CN);
        copyNumber = createCopyNumber(CN, correspondingBAF);
        variant = createVariant(alleleReadCount, VCN);

        biallelicProbability = SomaticPurityEnrichment.calculateBiallelic(copyNumber, variant);
        assertEquals(expectedBiallelicProbability, biallelicProbability, 0.01);

        /* Scenario 3: Opposing evidence -> expect probability of biallelic about 0.5 */

        // input values:
        CN = 2.06;
        MACN = 0.04;
        VCN = 1.57;
        alleleReadCount = 11;

        // expected output:
        expectedBiallelicProbability = 0.54;

        correspondingBAF = 1 - (MACN / CN);
        copyNumber = createCopyNumber(CN, correspondingBAF);
        variant = createVariant(alleleReadCount, VCN);

        biallelicProbability = SomaticPurityEnrichment.calculateBiallelic(copyNumber, variant);
        assertEquals(expectedBiallelicProbability, biallelicProbability, 0.01);

        /* Scenario 4 - another example of uncertain biallelic status */

        // input values:
        CN = 1.91;
        MACN = 0.0;
        VCN = 1.49;
        alleleReadCount = 16;

        // expected output:
        expectedBiallelicProbability = 0.533;

        correspondingBAF = 1 - (MACN / CN);
        copyNumber = createCopyNumber(CN, correspondingBAF);
        variant = createVariant(alleleReadCount, VCN);

        biallelicProbability = SomaticPurityEnrichment.calculateBiallelic(copyNumber, variant);
        assertEquals(expectedBiallelicProbability, biallelicProbability, 0.01);

        /* Scenario 5 - another example of uncertain biallelic status */

        // input values:
        CN = 2.16;
        MACN = 0.0;
        VCN = 1.51;
        alleleReadCount = 22;

        // expected output:
        expectedBiallelicProbability = 0.528;

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