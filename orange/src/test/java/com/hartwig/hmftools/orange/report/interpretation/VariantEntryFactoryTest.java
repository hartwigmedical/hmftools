package com.hartwig.hmftools.orange.report.interpretation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogTestFactory;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.orange.algo.purple.PurpleVariant;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleVariantFactory;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class VariantEntryFactoryTest {

    private static final double EPSILON = 1.0E-10;

    @Test
    public void canCreateVariantEntries() {
        PurpleVariant driverVariant = TestPurpleVariantFactory.builder()
                .reported(true)
                .gene("gene 1")
                .canonicalImpact(TestPurpleVariantFactory.impactBuilder().transcript("transcript 1").hgvsProteinImpact("impact 1").build())
                .addOtherImpacts(TestPurpleVariantFactory.impactBuilder().transcript("transcript 2").hgvsProteinImpact("impact 2").build())
                .adjustedCopyNumber(2)
                .adjustedVAF(1.3)
                .minorAlleleCopyNumber(1.2)
                .biallelic(false)
                .hotspot(Hotspot.NEAR_HOTSPOT)
                .subclonalLikelihood(0.3)
                .localPhaseSets(Lists.newArrayList(1))
                .build();

        PurpleVariant nonDriverVariant = TestPurpleVariantFactory.builder()
                .gene("gene 2")
                .canonicalImpact(TestPurpleVariantFactory.impactBuilder().hgvsProteinImpact("impact 3").build())
                .build();

        DriverCatalog canonicalDriver = DriverCatalogTestFactory.builder()
                .driver(DriverType.MUTATION)
                .gene("gene 1")
                .transcript("transcript 1")
                .isCanonical(true)
                .driverLikelihood(0.5)
                .build();

        DriverCatalog nonCanonicalDriver = DriverCatalogTestFactory.builder()
                .driver(DriverType.MUTATION)
                .gene("gene 1")
                .transcript("transcript 2")
                .isCanonical(false)
                .driverLikelihood(0.4)
                .build();

        List<PurpleVariant> variants = Lists.newArrayList(driverVariant, nonDriverVariant);
        List<DriverCatalog> drivers = Lists.newArrayList(canonicalDriver, nonCanonicalDriver);

        List<VariantEntry> entries = VariantEntryFactory.create(variants, drivers);

        assertEquals(3, entries.size());
        VariantEntry entry1 = findByGeneAndImpact(entries, "gene 1", "impact 1");
        assertTrue(entry1.isCanonical());
        assertEquals(2D, entry1.variantCopyNumber(), EPSILON);
        assertEquals(2D, entry1.totalCopyNumber(), EPSILON);
        assertEquals(1.2, entry1.minorAlleleCopyNumber(), EPSILON);
        assertFalse(entry1.biallelic());
        assertEquals(Hotspot.NEAR_HOTSPOT, entry1.hotspot());
        assertEquals(0.5, entry1.driverLikelihood(), EPSILON);
        assertEquals(0.7, entry1.clonalLikelihood(), EPSILON);

        VariantEntry entry2 = findByGeneAndImpact(entries, "gene 1", "impact 2");
        assertFalse(entry2.isCanonical());
        assertEquals(0.4, entry2.driverLikelihood(), EPSILON);

        VariantEntry entry3 = findByGeneAndImpact(entries, "gene 2", "impact 3");
        assertNull(entry3.driverLikelihood());
    }

    @NotNull
    private static VariantEntry findByGeneAndImpact(@NotNull List<VariantEntry> entries, @NotNull String geneToFind,
            @NotNull String impactToFind) {
        for (VariantEntry entry : entries) {
            if (entry.gene().equals(geneToFind) && entry.impact().equals(impactToFind)) {
                return entry;
            }
        }

        throw new IllegalStateException("Could not find variant entry with gene and " + geneToFind + " and impact " + impactToFind);
    }

    @Test
    public void canDetermineTranscriptImpact() {
        assertEquals("p.G12C",
                VariantEntryFactory.determineImpact(TestPurpleVariantFactory.impactBuilder()
                        .hgvsCodingImpact("c.123A>C")
                        .hgvsProteinImpact("p.Gly12Cys")
                        .addEffects(VariantEffect.MISSENSE)
                        .codingEffect(CodingEffect.MISSENSE)
                        .build()));
        assertEquals("c.123A>C splice",
                VariantEntryFactory.determineImpact(TestPurpleVariantFactory.impactBuilder()
                        .hgvsCodingImpact("c.123A>C")
                        .hgvsProteinImpact("p.?")
                        .addEffects(VariantEffect.MISSENSE)
                        .codingEffect(CodingEffect.SPLICE)
                        .build()));
        assertEquals("c.123A>C",
                VariantEntryFactory.determineImpact(TestPurpleVariantFactory.impactBuilder()
                        .hgvsCodingImpact("c.123A>C")
                        .hgvsProteinImpact(Strings.EMPTY)
                        .addEffects(VariantEffect.MISSENSE)
                        .codingEffect(CodingEffect.MISSENSE)
                        .build()));
        assertEquals("missense_variant",
                VariantEntryFactory.determineImpact(TestPurpleVariantFactory.impactBuilder()
                        .hgvsCodingImpact(Strings.EMPTY)
                        .hgvsProteinImpact(Strings.EMPTY)
                        .addEffects(VariantEffect.MISSENSE)
                        .codingEffect(CodingEffect.MISSENSE)
                        .build()));
    }
}