package com.hartwig.hmftools.orange.algo.util;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;

import org.junit.Test;

public class VariantTranscriptImpactCleanerTest
{
    @Test
    public void testCleanBrackets()
    {
        var impact = new VariantTranscriptImpact(
                "[ENSG00000109265", // this bracket needs to be removed
                "CRACD",
                "ENST00000264229",
                "synonymous_variant",
                false,
                "c.2187C>T",
                "p.Ser729=]"); // this bracket needs to be removed

        // act
        var cleanedImpact = VariantTranscriptImpactCleaner.cleanFields(impact);

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
        var impact = new VariantTranscriptImpact(
                " ENSG00000109265", // this whitespace needs to be removed
                "CRACD",
                "ENST00000264229 ", // this whitespace needs to be removed
                "synonymous_variant",
                false,
                "c.2187C>T",
                "p.Ser729=]"); // this bracket needs to be removed

        // act
        var cleanedImpact = VariantTranscriptImpactCleaner.cleanFields(impact);

        assertEquals("ENSG00000109265", cleanedImpact.GeneId);
        assertEquals("CRACD", cleanedImpact.GeneName);
        assertEquals("ENST00000264229", cleanedImpact.Transcript);
        assertEquals("synonymous_variant", cleanedImpact.Effects);
        assertFalse(cleanedImpact.SpliceRegion);
        assertEquals("c.2187C>T", cleanedImpact.HgvsCoding);
        assertEquals("p.Ser729=", cleanedImpact.HgvsProtein);
    }
}
