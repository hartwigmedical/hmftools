package com.hartwig.hmftools.orange.report.interpretation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.finding.Variants;
import com.hartwig.hmftools.datamodel.purple.HotspotType;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleVariantFactory;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.datamodel.TestVariantEntryFactory;
import com.hartwig.hmftools.datamodel.finding.SmallVariant;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class VariantsTest
{
    @Test
    public void canSortVariantEntries()
    {
        SmallVariant entry1 = TestVariantEntryFactory.builder().gene("gene A").driverLikelihood(1D).affectedCodon(600).build();
        SmallVariant entry2 = TestVariantEntryFactory.builder().gene("gene A").driverLikelihood(1D).affectedCodon(700).build();
        SmallVariant entry3 = TestVariantEntryFactory.builder().gene("gene B").driverLikelihood(1D).affectedCodon(600).build();
        SmallVariant entry4 = TestVariantEntryFactory.builder().gene("gene B").driverLikelihood(0.5).affectedCodon(600).build();
        SmallVariant entry5 = TestVariantEntryFactory.builder().gene("gene B").driverLikelihood(null).affectedCodon(600).build();
        SmallVariant entry6 = TestVariantEntryFactory.builder().gene("gene B").driverLikelihood(null).affectedCodon(null).build();

        List<SmallVariant> entries = Lists.newArrayList(entry5, entry6, entry1, entry3, entry2, entry4);
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
        SmallVariant canonical = TestVariantEntryFactory.builder().gene("gene").isCanonical(true).impact("impact").build();
        assertEquals("gene impact", Variants.variantField(canonical));

        SmallVariant nonCanonical = TestVariantEntryFactory.builder().gene("gene").isCanonical(false).impact("impact").build();
        assertEquals("gene (alt) impact", Variants.variantField(nonCanonical));
    }

    @Test
    public void canHandleAllHotspotValues()
    {
        for(HotspotType hotspot : HotspotType.values())
        {
            SmallVariant entry = TestVariantEntryFactory.builder().hotspot(hotspot).build();
            assertNotNull(Variants.hotspotField(entry));
        }
    }

    @Test
    public void canRenderBiallelicStatus()
    {
        SmallVariant biallelic = TestVariantEntryFactory.builder().biallelicProbability(1.).build();
        assertEquals("100%", Variants.biallelicLikelihoodField(biallelic));

        SmallVariant nonBiallelic = TestVariantEntryFactory.builder().biallelicProbability(0.2).build();
        assertEquals("20%", Variants.biallelicLikelihoodField(nonBiallelic));
    }

    @Test
    public void canRenderDriverLikelihood()
    {
        SmallVariant driver = TestVariantEntryFactory.builder().driverLikelihood(0.4).build();
        assertEquals("40%", Variants.driverLikelihoodField(driver));

        SmallVariant nonDriver = TestVariantEntryFactory.builder().driverLikelihood(null).build();
        assertEquals(Strings.EMPTY, Variants.driverLikelihoodField(nonDriver));
    }

    @Test
    public void canRenderClonalLikelihood()
    {
        SmallVariant clonal = TestVariantEntryFactory.builder().clonalLikelihood(1D).build();
        assertEquals("100%", Variants.clonalLikelihoodField(clonal));

        SmallVariant subclonal = TestVariantEntryFactory.builder().clonalLikelihood(0.3).build();
        assertEquals("30%", Variants.clonalLikelihoodField(subclonal));
    }

    @Test
    public void canRenderRnaDepthField()
    {
        SmallVariant missingRna = TestVariantEntryFactory.builder().rnaDepth(null).build();
        assertEquals(ReportResources.NOT_AVAILABLE, Variants.rnaDepthField(missingRna));

        SmallVariant proper = TestVariantEntryFactory.builder()
                .rnaDepth(TestPurpleVariantFactory.depthBuilder().alleleReadCount(10).totalReadCount(20).build())
                .build();
        assertEquals("10/20 (50%)", Variants.rnaDepthField(proper));

        SmallVariant noDepth = TestVariantEntryFactory.builder()
                .rnaDepth(TestPurpleVariantFactory.depthBuilder().alleleReadCount(0).totalReadCount(0).build())
                .build();
        assertEquals("0/0", Variants.rnaDepthField(noDepth));
    }

    @Test
    public void canRenderPhaseSetField()
    {
        SmallVariant missingPhaseSet = TestVariantEntryFactory.builder().localPhaseSets(null).build();
        assertEquals(Strings.EMPTY, Variants.phaseSetField(missingPhaseSet));

        SmallVariant emptyPhaseSet = TestVariantEntryFactory.builder().localPhaseSets(Lists.newArrayList()).build();
        assertEquals(Strings.EMPTY, Variants.phaseSetField(emptyPhaseSet));

        SmallVariant multiPhaseSet = TestVariantEntryFactory.builder().localPhaseSets(Lists.newArrayList(1, 2)).build();
        assertEquals("1, 2", Variants.phaseSetField(multiPhaseSet));
    }
}