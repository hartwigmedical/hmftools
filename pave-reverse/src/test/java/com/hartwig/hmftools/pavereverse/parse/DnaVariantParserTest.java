package com.hartwig.hmftools.pavereverse.parse;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.pavereverse.ReversePaveTestBase;
import com.hartwig.hmftools.pavereverse.dna.DnaVariant;

import org.junit.Before;
import org.junit.Test;

public class DnaVariantParserTest extends ReversePaveTestBase
{
    private DnaVariantParser parser;

    @Before
    public void setUp()
    {
        parser = reversePave.dnaVariantParser();
    }

    @Test
    public void parseSubstitutionTest()
    {
        final String transcriptId = "ENST00000322764";
        DnaVariant variant = parser.parse("ZYX", transcriptId, "6G>A");
        assertEquals(transcriptId, variant.Transcript.TransName);
        assertEquals("ZYX", variant.Gene.GeneName);
        assertEquals("G", variant.Ref);
        assertEquals("A", variant.Alt);
        assertEquals(6, variant.Position);
    }


}
