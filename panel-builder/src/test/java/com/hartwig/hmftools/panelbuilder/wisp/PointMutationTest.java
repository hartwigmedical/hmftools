package com.hartwig.hmftools.panelbuilder.wisp;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser.writeImpactDetails;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;

import org.junit.Test;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

// TODO

public class PointMutationTest
{
    @Test
    public void testPointMutationSequences()
    {
        TestUtils.MOCK_REF_GENOME.RefGenomeMap.put(CHR_1, TestUtils.REF_BASES_CHR_1);

        // an SNV
        int position = 20;
        String ref = TestUtils.REF_BASES_CHR_1.substring(position, position + 1);
        String alt = "A";
        SomaticMutation var = new SomaticMutation(createVariantContext(CHR_1, position, ref, alt), TestUtils.TEST_SAMPLE_ID);

        var.generateSequences(TestUtils.MOCK_REF_GENOME, TestUtils.TEST_CONFIG);

        String sequence = TestUtils.REF_BASES_CHR_1.substring(10, 20) + alt + TestUtils.REF_BASES_CHR_1.substring(21, 30);
        assertEquals(var.sequence(), sequence);

        // an MNV
        ref = TestUtils.REF_BASES_CHR_1.substring(position, position + 2);
        alt = "AA";
        var = new SomaticMutation(createVariantContext(CHR_1, position, ref, alt), TestUtils.TEST_SAMPLE_ID);

        var.generateSequences(TestUtils.MOCK_REF_GENOME, TestUtils.TEST_CONFIG);

        sequence = TestUtils.REF_BASES_CHR_1.substring(11, 20) + alt + TestUtils.REF_BASES_CHR_1.substring(22, 31);
        assertEquals(var.sequence(), sequence);

        // 3-base MNV
        ref = TestUtils.REF_BASES_CHR_1.substring(position, position + 3);
        alt = "AAA";
        var = new SomaticMutation(createVariantContext(CHR_1, position, ref, alt), TestUtils.TEST_SAMPLE_ID);

        var.generateSequences(TestUtils.MOCK_REF_GENOME, TestUtils.TEST_CONFIG);

        sequence = TestUtils.REF_BASES_CHR_1.substring(11, 20) + alt + TestUtils.REF_BASES_CHR_1.substring(23, 31);
        assertEquals(sequence, var.sequence());

        // 5-base delete
        int delLength = 5;
        ref = TestUtils.REF_BASES_CHR_1.substring(position, position + delLength);
        alt = ref.substring(0, 1);
        var = new SomaticMutation(createVariantContext(CHR_1, position, ref, alt), TestUtils.TEST_SAMPLE_ID);

        var.generateSequences(TestUtils.MOCK_REF_GENOME, TestUtils.TEST_CONFIG);

        sequence = TestUtils.REF_BASES_CHR_1.substring(10, 20) + alt + TestUtils.REF_BASES_CHR_1.substring(
                position + delLength, position + delLength + 9);
        assertEquals(sequence, var.sequence());

        // 10-base delete
        delLength = 10;
        ref = TestUtils.REF_BASES_CHR_1.substring(position, position + delLength);
        alt = ref.substring(0, 1);
        var = new SomaticMutation(createVariantContext(CHR_1, position, ref, alt), TestUtils.TEST_SAMPLE_ID);

        var.generateSequences(TestUtils.MOCK_REF_GENOME, TestUtils.TEST_CONFIG);

        sequence = TestUtils.REF_BASES_CHR_1.substring(10, 20) + alt + TestUtils.REF_BASES_CHR_1.substring(
                position + delLength, position + delLength + 9);
        assertEquals(sequence, var.sequence());

        // 2 base insert
        ref = TestUtils.REF_BASES_CHR_1.substring(position, position + 1);
        alt = ref + "AA";
        var = new SomaticMutation(createVariantContext(CHR_1, position, ref, alt), TestUtils.TEST_SAMPLE_ID);

        var.generateSequences(TestUtils.MOCK_REF_GENOME, TestUtils.TEST_CONFIG);

        sequence = TestUtils.REF_BASES_CHR_1.substring(11, 20) + alt + TestUtils.REF_BASES_CHR_1.substring(position + 1, position + 9);
        assertEquals(sequence, var.sequence());

        // 7 base insert
        ref = TestUtils.REF_BASES_CHR_1.substring(position, position + 1);
        alt = ref + "AAAAAAA";
        var = new SomaticMutation(createVariantContext(CHR_1, position, ref, alt), TestUtils.TEST_SAMPLE_ID);

        var.generateSequences(TestUtils.MOCK_REF_GENOME, TestUtils.TEST_CONFIG);

        sequence = TestUtils.REF_BASES_CHR_1.substring(14, 20) + alt + TestUtils.REF_BASES_CHR_1.substring(position + 1, position + 7);
        assertEquals(sequence, var.sequence());
    }

    private static VariantContext createVariantContext(final String chromosome, int position, final String ref, final String alt)
    {
        return createVariantContext(
                chromosome, position, ref, alt,
                new VariantImpact("", "", "", CodingEffect.NONE,
                        "", "", false, "",
                        CodingEffect.NONE, 0), TestUtils.DEFAULT_QUAL, 10);
    }

    private static VariantContext createVariantContext(
            final String chromosome, int position, final String ref, final String alt,
            final VariantImpact variantImpact, double qual, int tumorFragments)
    {
        VariantContextBuilder builder = new VariantContextBuilder();

        List<Allele> alleles = Lists.newArrayList();

        alleles.add(Allele.create(ref, true));
        alleles.add(Allele.create(alt, false));

        Map<String, Object> refAttributes = Maps.newHashMap();
        Map<String, Object> tumorAttributes = Maps.newHashMap();

        Map<String, Object> commonAttributes = Maps.newHashMap();

        Genotype gtNormal = new GenotypeBuilder()
                .attributes(refAttributes)
                .name(TestUtils.TEST_REF_ID)
                .DP(-1)
                .noAD()
                .noPL()
                .GQ(-1)
                .make();

        Genotype gtTumor = new GenotypeBuilder()
                .attributes(tumorAttributes)
                .name(TestUtils.TEST_SAMPLE_ID)
                .DP(-1)
                .AD(new int[] { 0, tumorFragments })
                .noPL()
                .GQ(-1)
                .make();

        GenotypesContext genotypesContext = GenotypesContext.create(gtNormal, gtTumor);

        // String filters = PASS;

        double logError = -(qual / 10.0);

        // int positionEnd =

        VariantContext variantContext = builder
                .source("SOURCE")
                .chr(chromosome)
                .start(position)
                .stop(position + ref.length() - 1)
                .alleles(alleles)
                .genotypes(genotypesContext)
                .attributes(commonAttributes)
                .log10PError(logError)
                .unfiltered()
                .make(true);

        writeImpactDetails(variantContext, variantImpact);

        return variantContext;
    }
}
