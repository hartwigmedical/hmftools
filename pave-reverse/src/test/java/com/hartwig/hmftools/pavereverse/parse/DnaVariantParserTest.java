package com.hartwig.hmftools.pavereverse.parse;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.pavereverse.ReversePaveTestBase;
import com.hartwig.hmftools.pavereverse.dna.DnaVariant;
import com.hartwig.hmftools.pavereverse.dna.InExonDownstreamOfCodingEnd;
import com.hartwig.hmftools.pavereverse.dna.InExon;
import com.hartwig.hmftools.pavereverse.dna.InIntronAfterExon;
import com.hartwig.hmftools.pavereverse.dna.InIntronBeforeExon;
import com.hartwig.hmftools.pavereverse.dna.InExonUpstreamOfCodingStart;
import com.hartwig.hmftools.pavereverse.dna.InIntronDownstreamOfCodingEnd;
import com.hartwig.hmftools.pavereverse.dna.InIntronUpstreamOfCodingStart;

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
    public void parseSubstitutionInExon()
    {
        DnaVariant variant = parser.parse(zyx, zyxCanonical, "6G>A");
        assertEquals(zyxCanonical, variant.transcriptName());
        assertEquals("ZYX", variant.geneName());
        assertEquals("G", variant.Ref);
        assertEquals("A", variant.Alt);
        assertEquals(6, ((InExon) variant.Address).IndexOfBaseInCodingBases);
    }

    @Test
    public void parseSubstitutionUpstreamOfCodingStart()
    {
        DnaVariant variant = parser.parse("ZYX", zyxCanonical, "c.-6G>A");
        assertEquals(zyxCanonical, variant.transcriptName());
        assertEquals(zyx, variant.geneName());
        assertEquals("G", variant.Ref);
        assertEquals("A", variant.Alt);
        assertEquals(-6, ((InExonUpstreamOfCodingStart) variant.Address).IndexUpstreamOfStart);
    }

    @Test
    public void parseSubstitutionDownstreamOfCodingEnd()
    {
        DnaVariant variant = parser.parse(zyx, zyxCanonical, "*12C>G");
        assertEquals(zyxCanonical, variant.transcriptName());
        assertEquals(zyx, variant.geneName());
        assertEquals("C", variant.Ref);
        assertEquals("G", variant.Alt);
        assertEquals(12, ((InExonDownstreamOfCodingEnd) variant.Address).IndexDownstreamOfEnd);
    }

    @Test
    public void parseSubstitutionInIntron()
    {
        DnaVariant variant = parser.parse(zyx, zyxCanonical, "c.208+1G>A");
        assertEquals("G", variant.Ref);
        assertEquals("A", variant.Alt);
        assertEquals(208, ((InIntronAfterExon) variant.Address).ExonBaseIndex);
        assertEquals(1, ((InIntronAfterExon) variant.Address).IndexAfterExonBase);

        variant = parser.parse(zyx, zyxCanonical, "209-1G>A");
        assertEquals("G", variant.Ref);
        assertEquals("A", variant.Alt);
        assertEquals(209, ((InIntronBeforeExon) variant.Address).ExonBaseIndex);
        assertEquals(1, ((InIntronBeforeExon) variant.Address).IndexBeforeExonBase);
    }

    @Test
    public void parseSubstitutionInIntronUpstreamOfCodingStart()
    {
        DnaVariant variant = parser.parse(zyx, zyxCanonical, "c.-16-1G>A");
        assertEquals("G", variant.Ref);
        assertEquals("A", variant.Alt);
        assertEquals(-16, ((InIntronUpstreamOfCodingStart) variant.Address).IndexOfExonicBase);
        assertEquals(-1, ((InIntronUpstreamOfCodingStart) variant.Address).RelativePositionOfIntronicBase);

        variant = parser.parse(zyx, zyxCanonical, "-17+1G>A");
        assertEquals("G", variant.Ref);
        assertEquals("A", variant.Alt);
        assertEquals(-17, ((InIntronUpstreamOfCodingStart) variant.Address).IndexOfExonicBase);
        assertEquals(1, ((InIntronUpstreamOfCodingStart) variant.Address).RelativePositionOfIntronicBase);
    }

    @Test
    public void parseSubstitutionInIntronDownstreamOfCodingEnd()
    {
        DnaVariant variant = parser.parse(tatdn2, tatdn2Canonical, "c.*38+1G>A");
        assertEquals("G", variant.Ref);
        assertEquals("A", variant.Alt);
        assertEquals(38, ((InIntronDownstreamOfCodingEnd) variant.Address).IndexOfExonicBase);
        assertEquals(1, ((InIntronDownstreamOfCodingEnd) variant.Address).RelativePositionOfIntronicBase);

        variant = parser.parse(tatdn2, tatdn2Canonical, "*39-3C>A");
        assertEquals("C", variant.Ref);
        assertEquals("A", variant.Alt);
        assertEquals(39, ((InIntronDownstreamOfCodingEnd) variant.Address).IndexOfExonicBase);
        assertEquals(-3, ((InIntronDownstreamOfCodingEnd) variant.Address).RelativePositionOfIntronicBase);
    }
}
