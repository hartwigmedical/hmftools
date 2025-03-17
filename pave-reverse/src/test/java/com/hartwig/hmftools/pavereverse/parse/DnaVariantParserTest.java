package com.hartwig.hmftools.pavereverse.parse;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.pavereverse.ReversePaveTestBase;
import com.hartwig.hmftools.pavereverse.dna.DnaVariant;
import com.hartwig.hmftools.pavereverse.dna.DownstreamOfCodingEndAddress;
import com.hartwig.hmftools.pavereverse.dna.InExonAddress;
import com.hartwig.hmftools.pavereverse.dna.UpstreamOfCodingStartAddress;

import org.junit.Before;
import org.junit.Test;

public class DnaVariantParserTest extends ReversePaveTestBase
{
    private DnaVariantParser parser;
    private final String zyx = "ZYX";
    private final String zyxCanonical = "ENST00000322764";

    @Before
    public void setUp()
    {
        parser = reversePave.dnaVariantParser();
    }

    @Test
    public void parseSubstitutionInExonTest()
    {
        DnaVariant variant = parser.parse(zyx, zyxCanonical, "6G>A");
        assertEquals(zyxCanonical, variant.transcriptName());
        assertEquals("ZYX", variant.geneName());
        assertEquals("G", variant.Ref);
        assertEquals("A", variant.Alt);
        assertEquals(6, ((InExonAddress) variant.Address).IndexOfBaseInCodingBases);
    }

    @Test
    public void parseSubstitutionUpstreamOfCodingStartTest()
    {
        DnaVariant variant = parser.parse("ZYX", zyxCanonical, "-6G>A");
        assertEquals(zyxCanonical, variant.transcriptName());
        assertEquals(zyx, variant.geneName());
        assertEquals("G", variant.Ref);
        assertEquals("A", variant.Alt);
        assertEquals(6, ((UpstreamOfCodingStartAddress) variant.Address).IndexUpstreamOfStart);
    }

    @Test
    public void parseDownstreamOfCodingEndTest()
    {
        DnaVariant variant = parser.parse(zyx, zyxCanonical, "*12C>G");
        assertEquals(zyxCanonical, variant.transcriptName());
        assertEquals(zyx, variant.geneName());
        assertEquals("C", variant.Ref);
        assertEquals("G", variant.Alt);
        assertEquals(12, ((DownstreamOfCodingEndAddress) variant.Address).IndexDownstreamOfEnd);
    }
}
