package com.hartwig.hmftools.orange.algo.util;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;

import org.junit.Test;

public class VariantTranscriptImpactCleanerTest
{
    @Test
    public void testCleanBrackets()
    {
        VariantTranscriptImpact impact = new VariantTranscriptImpact(
                "[ENSG00000109265", // this bracket needs to be removed
                "CRACD",
                "ENST00000264229",
                "synonymous_variant",
                false,
                "c.2187C>T",
                "p.Ser729=]", "", null, null); // this bracket needs to be removed

        VariantTranscriptImpact cleanedImpact = VariantTranscriptImpactCleaner.cleanFields(impact);

        assertEquals("ENSG00000109265", cleanedImpact.GeneId);
        assertEquals("CRACD", cleanedImpact.GeneName);
        assertEquals("ENST00000264229", cleanedImpact.Transcript);
        assertEquals("synonymous_variant", cleanedImpact.Effects);
        assertFalse(cleanedImpact.SpliceRegion);
        assertEquals("c.2187C>T", cleanedImpact.HgvsCoding);
        assertEquals("p.Ser729=", cleanedImpact.HgvsProtein);
    }

    @Test
    public void testTrailingSpacesAndBracket()
    {
        VariantTranscriptImpact impact = new VariantTranscriptImpact(
                " ENSG00000109265", // this whitespace needs to be removed
                "CRACD",
                "ENST00000264229 ", // this whitespace needs to be removed
                "synonymous_variant",
                false,
                "c.2187C>T",
                "p.Ser729=]", "", null, null); // this bracket needs to be removed

        VariantTranscriptImpact cleanedImpact = VariantTranscriptImpactCleaner.cleanFields(impact);

        assertEquals("ENSG00000109265", cleanedImpact.GeneId);
        assertEquals("CRACD", cleanedImpact.GeneName);
        assertEquals("ENST00000264229", cleanedImpact.Transcript);
        assertEquals("synonymous_variant", cleanedImpact.Effects);
        assertFalse(cleanedImpact.SpliceRegion);
        assertEquals("c.2187C>T", cleanedImpact.HgvsCoding);
        assertEquals("p.Ser729=", cleanedImpact.HgvsProtein);
    }

    @Test
    public void testImpactWithNullFields()
    {
        VariantTranscriptImpact impact = new VariantTranscriptImpact(
                null,
                null,
                null,
                null,
                true,
                null,
                null, "", null, null);

        VariantTranscriptImpact cleanedImpact = VariantTranscriptImpactCleaner.cleanFields(impact);

        assertNull(cleanedImpact.GeneId);
        assertNull(cleanedImpact.GeneName);
        assertNull(cleanedImpact.Transcript);
        assertNull(cleanedImpact.Effects);
        assertTrue(cleanedImpact.SpliceRegion);
        assertNull(cleanedImpact.HgvsCoding);
        assertNull(cleanedImpact.HgvsProtein);
    }
}
