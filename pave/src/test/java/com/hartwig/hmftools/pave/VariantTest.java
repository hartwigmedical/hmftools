package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_VARIANT_CN;
import static com.hartwig.hmftools.common.variant.SageVcfTags.MICROHOMOLOGY;
import static com.hartwig.hmftools.common.variant.SageVcfTags.REPEAT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.REPEAT_SEQUENCE;
import static com.hartwig.hmftools.common.variant.VariantTier.HOTSPOT;
import static com.hartwig.hmftools.common.variant.VariantTier.PANEL;
import static com.hartwig.hmftools.common.variant.VariantTier.TIER;
import static com.hartwig.hmftools.pave.ChromosomeTask.applyFilters;
import static com.hartwig.hmftools.pave.annotation.GnomadAnnotation.PON_GNOMAD_FILTER;
import static com.hartwig.hmftools.pave.annotation.PonAnnotation.PON_ARTEFACT_FILTER;
import static com.hartwig.hmftools.pave.annotation.PonAnnotation.PON_FILTER;

import static htsjdk.variant.vcf.VCFConstants.ALLELE_FREQUENCY_KEY;
import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.StringCache;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.pave.annotation.PonAnnotation;
import com.hartwig.hmftools.pave.annotation.PonChrCache;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class VariantTest
{
    @Test
    public void testPointMutation()
    {
        // SNVs and MNVs
        VariantData var = new VariantData(CHR_1, 100, "A", "G");

        assertTrue(var.isBaseChange());
        assertEquals(100, var.EndPosition);
        assertEquals(1, var.altPositions().size());
        assertTrue(var.altBasesAbove(99));
        assertFalse(var.altBasesAbove(100));
        assertTrue(var.altBasesBelow(101));
        assertFalse(var.altBasesBelow(100));

        assertTrue(var.altPositionsOverlap(90, 100));
        assertFalse(var.altPositionsOverlap(90, 99));
        assertTrue(var.altPositionsOverlap(100, 110));
        assertFalse(var.altPositionsOverlap(101, 110));

        assertTrue(var.altPositionsWithin(100, 110));
        assertTrue(var.altPositionsWithin(90, 100));

        var = new VariantData(CHR_1, 100, "AGA", "GCG");

        assertTrue(var.isBaseChange());
        assertEquals(102, var.EndPosition);
        assertEquals(3, var.altPositions().size());
        assertTrue(var.altBasesAbove(99));
        assertFalse(var.altBasesAbove(100));
        assertFalse(var.altBasesAbove(102));
        assertTrue(var.altBasesBelow(103));
        assertFalse(var.altBasesBelow(102));
        assertFalse(var.altBasesBelow(100));

        assertTrue(var.altPositionsOverlap(90, 100));
        assertFalse(var.altPositionsOverlap(103, 110));

        assertTrue(var.altPositionsWithin(100, 110));
        assertFalse(var.altPositionsWithin(102, 110));
        assertTrue(var.altPositionsWithin(90, 102));
        assertFalse(var.altPositionsWithin(90, 101));
    }

    @Test
    public void testDeletion()
    {
        VariantData var = new VariantData(CHR_1, 100, "ACGT", "A");

        assertTrue(var.isDeletion());
        assertEquals(104, var.EndPosition);
        assertEquals(3, var.altPositions().size());
        assertTrue(var.altBasesAbove(99));
        assertTrue(var.altBasesAbove(100));
        assertFalse(var.altBasesAbove(101));
        assertTrue(var.altBasesBelow(104));
        assertFalse(var.altBasesBelow(103));

        assertFalse(var.altPositionsOverlap(90, 100));
        assertTrue(var.altPositionsOverlap(101, 110));
        assertFalse(var.altPositionsOverlap(104, 110));
        assertTrue(var.altPositionsOverlap(103, 110));

        assertTrue(var.altPositionsWithin(100, 110));
        assertTrue(var.altPositionsWithin(101, 110));
        assertFalse(var.altPositionsWithin(102, 110));
        assertTrue(var.altPositionsWithin(90, 104));
        assertTrue(var.altPositionsWithin(90, 103));
        assertFalse(var.altPositionsWithin(90, 102));
    }

    @Test
    public void testInsertion()
    {
        VariantData var = new VariantData(CHR_1, 100, "A", "AAAAA"); // 4 bases

        assertTrue(var.isInsert());
        assertEquals(101, var.EndPosition);
        assertTrue(var.altPositions().isEmpty());

        assertTrue(var.altBasesAbove(100));
        assertFalse(var.altBasesAbove(101));
        assertTrue(var.altBasesBelow(101));
        assertFalse(var.altBasesBelow(100));

        assertFalse(var.altPositionsOverlap(90, 100));
        assertTrue(var.altPositionsOverlap(90, 101));
        assertFalse(var.altPositionsOverlap(101, 110));
        assertTrue(var.altPositionsOverlap(100, 110));

        assertTrue(var.altPositionsWithin(100, 110));
        assertFalse(var.altPositionsWithin(101, 110));
        assertTrue(var.altPositionsWithin(90, 101));
        assertFalse(var.altPositionsWithin(90, 100));
    }

    private static final String SAMPLE_ID = "SAMPLE";

    @Test
    public void testPonFilters()
    {
        VariantData var = new VariantData(CHR_1, 100, "A", "T");

        // first test a hotspot variant
        VariantContext variantContext = buildVariantContext(var, HOTSPOT, 0.2, 3, "");
        var.setContext(variantContext);

        PonAnnotation standardPon = new PonAnnotation(null, false);
        PonChrCache artefactsPon = new PonChrCache(CHR_1, new StringCache());

        applyFilters(var, SAMPLE_ID, standardPon, artefactsPon);
        assertFilters(var, false, false, false);

        // Gnomad
        var.setGnomadFrequency(0.02);
        applyFilters(var, SAMPLE_ID, standardPon, artefactsPon);
        assertFilters(var, false, false, true);

        // standard PON
        var.setPonFrequency(11, 7, 100);
        applyFilters(var, SAMPLE_ID, standardPon, artefactsPon);
        assertFilters(var, true, false, true);

        // artefact PON
        artefactsPon.addEntry(var.Position, var.Ref, var.Alt, 11, 11, 100);
        applyFilters(var, SAMPLE_ID, standardPon, artefactsPon);
        assertFilters(var, true, true, true);

        artefactsPon.clear();
        artefactsPon.addEntry(var.Position, var.Ref, var.Alt, 5, 11, 20);
        applyFilters(var, SAMPLE_ID, standardPon, artefactsPon);
        assertFilters(var, true, false, true);

        // non-hotspot variant
        variantContext = buildVariantContext(var, PANEL, 0.2, 3, "");
        var.setContext(variantContext);

        var.setGnomadFrequency(0.001);
        applyFilters(var, SAMPLE_ID, standardPon, artefactsPon);
        assertFilters(var, true, true, true);
    }

    private static void assertFilters(final VariantData var, boolean ponStandard, boolean ponArtefact, boolean gnomad)
    {
        if(ponStandard)
            assertTrue(var.filters().contains(PON_FILTER));
        else
            assertFalse(var.filters().contains(PON_FILTER));

        if(ponArtefact)
            assertTrue(var.filters().contains(PON_ARTEFACT_FILTER));
        else
            assertFalse(var.filters().contains(PON_ARTEFACT_FILTER));

        if(gnomad)
            assertTrue(var.filters().contains(PON_GNOMAD_FILTER));
        else
            assertFalse(var.filters().contains(PON_GNOMAD_FILTER));
    }

    private static VariantContext buildVariantContext(
            final VariantData variant, final VariantTier tier, double af, int repeatCount, final String repeatBases)
    {
        GenotypeBuilder genotypeBuilder = new GenotypeBuilder(SAMPLE_ID);
        genotypeBuilder.AD(new int[] { 10, 10 });
        genotypeBuilder.DP(100);
        genotypeBuilder.alleles(Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL));
        genotypeBuilder.attribute(ALLELE_FREQUENCY_KEY, 0.1);
        Genotype genotype = genotypeBuilder.make();

        List<Allele> alleles = List.of(Allele.create(variant.Ref, true), Allele.create(variant.Alt, false));

        VariantContextBuilder builder = new VariantContextBuilder();
        builder.chr(variant.Chromosome).start(variant.Position);
        builder.alleles(alleles);
        builder.computeEndFromAlleles(alleles, variant.Position);

        builder.genotypes(List.of(genotype));

        builder.attribute(TIER, tier);
        builder.attribute(REPEAT_SEQUENCE, repeatBases);
        builder.attribute(REPEAT_COUNT, repeatCount);

        return builder.make();
    }
}
