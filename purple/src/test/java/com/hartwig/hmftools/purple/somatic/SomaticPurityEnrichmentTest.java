package com.hartwig.hmftools.purple.somatic;


import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_VARIANT_CN;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_QUALITY;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNTS;
import static com.hartwig.hmftools.purple.TestUtils.REF_SAMPLE_ID;
import static com.hartwig.hmftools.purple.TestUtils.SAMPLE_ID;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.SegmentSupport;

import org.junit.Test;
import static junit.framework.TestCase.assertEquals;

import java.util.List;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import junit.framework.TestCase;

public class SomaticPurityEnrichmentTest extends TestCase
{
    @Test
    public void testCalculateBiallelic()
    {
        /* 
        input values:
            CN = 2.00
            MACN = 0.975
            VCN = 2.44
            AF = 1.2238
        
        Expected biallelic probability = 0.033
         */
        
        // MACN = (1 - BAF) * CN <=> BAF = 1 - MACN / CN
        double correspondingBAF = 1 - (0.975 / 2.00);
                
        PurpleCopyNumber copyNumber = createCopyNumber(2.00, correspondingBAF);

        SomaticVariant variant = createVariant(100, 10, 2);

        // assertEquals( 0.975, copyNumber.minorAlleleCopyNumber(), 0.01);
        // assertEquals(2.00, copyNumber.averageTumorCopyNumber(), 0.01);

        double biallelicValue = SomaticPurityEnrichment.calculateBiallelic(copyNumber, variant);
        // assertEquals(0.5, biallelicValue, 0.01);
    }

    private static SomaticVariant createVariant(int depth, int alleleSupport, double variantCopyNumber)
    {
        GenotypeBuilder genotypeBuilder = new GenotypeBuilder(SAMPLE_ID);
        genotypeBuilder.DP(depth).AD(new int[] { alleleSupport, alleleSupport });
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