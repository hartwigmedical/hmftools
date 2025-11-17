package com.hartwig.hmftools.orange.report.interpretation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.finding.Variants;
import com.hartwig.hmftools.datamodel.purple.HotspotType;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleVariantFactory;
import com.hartwig.hmftools.orange.algo.util.PurpleDriverTestFactory;
import com.hartwig.hmftools.datamodel.finding.SmallVariant;
import com.hartwig.hmftools.orange.report.finding.TestVariantEntryFactory;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class VariantsTest
{
    @Test
    public void canSortVariantEntries()
    {
        SmallVariant entry1 = TestVariantEntryFactory.builder("gene A")
                .driver(PurpleDriverTestFactory.builder().driverLikelihood(1D).build())
                .transcriptImpact(TestPurpleVariantFactory.impactBuilder().affectedCodon(600).build()).build();
        SmallVariant entry2 = TestVariantEntryFactory.builder("gene A")
                .driver(PurpleDriverTestFactory.builder().driverLikelihood(1D).build())
                .transcriptImpact(TestPurpleVariantFactory.impactBuilder().affectedCodon(700).build()).build();
        SmallVariant entry3 = TestVariantEntryFactory.builder("gene B")
                .driver(PurpleDriverTestFactory.builder().driverLikelihood(1D).build())
                .transcriptImpact(TestPurpleVariantFactory.impactBuilder().affectedCodon(600).build()).build();
        SmallVariant entry4 = TestVariantEntryFactory.builder("gene B")
                .driver(PurpleDriverTestFactory.builder().driverLikelihood(0.5).build())
                .transcriptImpact(TestPurpleVariantFactory.impactBuilder().affectedCodon(600).build()).build();
        SmallVariant entry5 = TestVariantEntryFactory.builder("gene B")
                .transcriptImpact(TestPurpleVariantFactory.impactBuilder().affectedCodon(600).build()).build();
        SmallVariant entry6 = TestVariantEntryFactory.builder("gene B")
                .transcriptImpact(TestPurpleVariantFactory.impactBuilder().affectedCodon(null).build()).build();

        List<SmallVariant> entries = List.of(entry5, entry6, entry1, entry3, entry2, entry4);
        List<SmallVariant> sorted = Variants.sort(entries);

        assertEquals(entry1, sorted.get(0));
        assertEquals(entry2, sorted.get(1));
        assertEquals(entry3, sorted.get(2));
        assertEquals(entry4, sorted.get(3));
        assertEquals(entry5, sorted.get(4));
        assertEquals(entry6, sorted.get(5));
    }

    @Test
    public void canRenderVariantField()
    {
        SmallVariant canonical = TestVariantEntryFactory.builder("gene").isCanonical(true)
                .transcriptImpact(TestPurpleVariantFactory.impactBuilder().hgvsProteinImpact("impact").build()).build();
        assertEquals("gene impact", Variants.variantField(canonical, true));

        SmallVariant nonCanonical = TestVariantEntryFactory.builder("gene").isCanonical(false)
                .transcriptImpact(TestPurpleVariantFactory.impactBuilder().hgvsProteinImpact("impact").build()).build();
        assertEquals("gene (alt) impact", Variants.variantField(nonCanonical, true));
    }

    @Test
    public void canHandleAllHotspotValues()
    {
        for(HotspotType hotspot : HotspotType.values())
        {
            SmallVariant entry = TestVariantEntryFactory.builder("gene")
                    .purpleVariant(TestPurpleVariantFactory.builder().hotspot(hotspot).build()).build();
            assertNotNull(Variants.hotspotField(entry));
        }
    }

    @Test
    public void canRenderBiallelicStatus()
    {
        SmallVariant biallelic = TestVariantEntryFactory.builder("gene")
                .purpleVariant(TestPurpleVariantFactory.builder().biallelicProbability(1.).build()).build();
        assertEquals("100%", Variants.biallelicLikelihoodField(biallelic));

        SmallVariant nonBiallelic = TestVariantEntryFactory.builder("gene")
                .purpleVariant(TestPurpleVariantFactory.builder().biallelicProbability(0.2).build()).build();
        assertEquals("20%", Variants.biallelicLikelihoodField(nonBiallelic));
    }

    @Test
    public void canRenderDriverLikelihood()
    {
        SmallVariant driver = TestVariantEntryFactory.builder("gene")
                .driver(PurpleDriverTestFactory.builder().driverLikelihood(0.4).build()).build();
        assertEquals("40%", Variants.driverLikelihoodField(driver));

        SmallVariant nonDriver = TestVariantEntryFactory.builder("gene")
                .driver(null).build();
        assertEquals(Strings.EMPTY, Variants.driverLikelihoodField(nonDriver));
    }

    @Test
    public void canRenderClonalLikelihood()
    {
        SmallVariant clonal = TestVariantEntryFactory.builder("gene")
                .purpleVariant(TestPurpleVariantFactory.builder().subclonalLikelihood(0.0).build()).build();
        assertEquals("100%", Variants.clonalLikelihoodField(clonal));

        SmallVariant subclonal = TestVariantEntryFactory.builder("gene")
                .purpleVariant(TestPurpleVariantFactory.builder().subclonalLikelihood(0.7).build()).build();
        assertEquals("30%", Variants.clonalLikelihoodField(subclonal));
    }

    @Test
    public void canRenderRnaDepthField()
    {
        String notAvailableString = "N/A";

        SmallVariant missingRna = TestVariantEntryFactory.builder("gene")
                .purpleVariant(TestPurpleVariantFactory.builder().rnaDepth(null).build()).build();
        assertEquals(notAvailableString, Variants.rnaDepthField(missingRna, notAvailableString));

        SmallVariant proper = TestVariantEntryFactory.builder("gene")
                .purpleVariant(TestPurpleVariantFactory.builder()
                .rnaDepth(TestPurpleVariantFactory.depthBuilder().alleleReadCount(10).totalReadCount(20).build())
                .build()).build();
        assertEquals("10/20 (50%)", Variants.rnaDepthField(proper, notAvailableString));

        SmallVariant noDepth = TestVariantEntryFactory.builder("gene")
                .purpleVariant(TestPurpleVariantFactory.builder()
                .rnaDepth(TestPurpleVariantFactory.depthBuilder().alleleReadCount(0).totalReadCount(0).build())
                .build()).build();
        assertEquals("0/0", Variants.rnaDepthField(noDepth, notAvailableString));
    }

    @Test
    public void canRenderPhaseSetField()
    {
        SmallVariant missingPhaseSet = TestVariantEntryFactory.builder("gene")
                .purpleVariant(TestPurpleVariantFactory.builder().localPhaseSets(null).build()).build();
        assertEquals(Strings.EMPTY, Variants.phaseSetField(missingPhaseSet));

        SmallVariant emptyPhaseSet = TestVariantEntryFactory.builder("gene")
                .purpleVariant(TestPurpleVariantFactory.builder().localPhaseSets(Lists.newArrayList()).build()).build();
        assertEquals(Strings.EMPTY, Variants.phaseSetField(emptyPhaseSet));

        SmallVariant multiPhaseSet = TestVariantEntryFactory.builder("gene")
                .purpleVariant(TestPurpleVariantFactory.builder().localPhaseSets(Lists.newArrayList(1, 2)).build()).build();
        assertEquals("1, 2", Variants.phaseSetField(multiPhaseSet));
    }
}