package com.hartwig.hmftools.pavereverse.parse;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.pavereverse.ReversePaveTestBase;
import com.hartwig.hmftools.pavereverse.dna.DeletionVariant;
import com.hartwig.hmftools.pavereverse.dna.DuplicationVariant;
import com.hartwig.hmftools.pavereverse.dna.InExonDownstreamOfCodingEnd;
import com.hartwig.hmftools.pavereverse.dna.InExon;
import com.hartwig.hmftools.pavereverse.dna.InIntronAfterExon;
import com.hartwig.hmftools.pavereverse.dna.InIntronBeforeExon;
import com.hartwig.hmftools.pavereverse.dna.InExonUpstreamOfCodingStart;
import com.hartwig.hmftools.pavereverse.dna.InIntronDownstreamOfCodingEnd;
import com.hartwig.hmftools.pavereverse.dna.InIntronUpstreamOfCodingStart;
import com.hartwig.hmftools.pavereverse.dna.InsertionVariant;
import com.hartwig.hmftools.pavereverse.dna.SubstitutionVariant;

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
    public void parseSubstitutionWithDeletionAndInsertionGiven()
    {
        SubstitutionVariant sub = (SubstitutionVariant) parser.parse(zyx, zyxCanonical, "c.3_7delGCGGCinsAA");
        assertEquals(3, ((InExon)sub.AddressOfChangeStart).IndexOfBaseInCodingBases);
        assertEquals(7, ((InExon)sub.AddressOfChangeEnd).IndexOfBaseInCodingBases);
        assertEquals("GCGGC", sub.Ref);
        assertEquals("AA", sub.Alt);
    }

    @Test
    public void parseComplexSubstitution()
    {
        SubstitutionVariant sub = (SubstitutionVariant) parser.parse(zyx, zyxCanonical, "c.3_10delinsAA");
        assertEquals(3, ((InExon)sub.AddressOfChangeStart).IndexOfBaseInCodingBases);
        assertEquals(10, ((InExon)sub.AddressOfChangeEnd).IndexOfBaseInCodingBases);
        assertEquals("", sub.Ref);
        assertEquals("AA", sub.Alt);
    }

    @Test
    public void parseSubstitutionInExon()
    {
        SubstitutionVariant variant = (SubstitutionVariant) parser.parse(zyx, zyxCanonical, "6G>A");
        assertEquals(zyxCanonical, variant.transcriptName());
        assertEquals("ZYX", variant.geneName());
        assertEquals("G", variant.Ref);
        assertEquals("A", variant.Alt);
        assertEquals(6, ((InExon) variant.AddressOfChangeStart).IndexOfBaseInCodingBases);
    }

    @Test
    public void parseSubstitutionUpstreamOfCodingStart()
    {
        SubstitutionVariant variant = (SubstitutionVariant) parser.parse("ZYX", zyxCanonical, "c.-6G>A");
        assertEquals(zyxCanonical, variant.transcriptName());
        assertEquals(zyx, variant.geneName());
        assertEquals("G", variant.Ref);
        assertEquals("A", variant.Alt);
        assertEquals(-6, ((InExonUpstreamOfCodingStart) variant.AddressOfChangeStart).IndexUpstreamOfStart);
    }

    @Test
    public void parseSubstitutionDownstreamOfCodingEnd()
    {
        SubstitutionVariant variant = (SubstitutionVariant) parser.parse(zyx, zyxCanonical, "*12C>G");
        assertEquals(zyxCanonical, variant.transcriptName());
        assertEquals(zyx, variant.geneName());
        assertEquals("C", variant.Ref);
        assertEquals("G", variant.Alt);
        assertEquals(12, ((InExonDownstreamOfCodingEnd) variant.AddressOfChangeStart).IndexDownstreamOfEnd);
    }

    @Test
    public void parseSubstitutionInIntron()
    {
        SubstitutionVariant variant = (SubstitutionVariant) parser.parse(zyx, zyxCanonical, "c.208+1G>A");
        assertEquals("G", variant.Ref);
        assertEquals("A", variant.Alt);
        assertEquals(208, ((InIntronAfterExon) variant.AddressOfChangeStart).ExonBaseIndex);
        assertEquals(1, ((InIntronAfterExon) variant.AddressOfChangeStart).IndexAfterExonBase);

        variant = (SubstitutionVariant) parser.parse(zyx, zyxCanonical, "209-1G>A");
        assertEquals("G", variant.Ref);
        assertEquals("A", variant.Alt);
        assertEquals(209, ((InIntronBeforeExon) variant.AddressOfChangeStart).ExonBaseIndex);
        assertEquals(1, ((InIntronBeforeExon) variant.AddressOfChangeStart).IndexBeforeExonBase);
    }

    @Test
    public void parseSubstitutionInIntronUpstreamOfCodingStart()
    {
        SubstitutionVariant variant = (SubstitutionVariant) parser.parse(zyx, zyxCanonical, "c.-16-1G>A");
        assertEquals("G", variant.Ref);
        assertEquals("A", variant.Alt);
        assertEquals(-16, ((InIntronUpstreamOfCodingStart) variant.AddressOfChangeStart).IndexOfExonicBase);
        assertEquals(-1, ((InIntronUpstreamOfCodingStart) variant.AddressOfChangeStart).RelativePositionOfIntronicBase);

        variant = (SubstitutionVariant) parser.parse(zyx, zyxCanonical, "-17+1G>A");
        assertEquals("G", variant.Ref);
        assertEquals("A", variant.Alt);
        assertEquals(-17, ((InIntronUpstreamOfCodingStart) variant.AddressOfChangeStart).IndexOfExonicBase);
        assertEquals(1, ((InIntronUpstreamOfCodingStart) variant.AddressOfChangeStart).RelativePositionOfIntronicBase);
    }

    @Test
    public void parseSubstitutionInIntronDownstreamOfCodingEnd()
    {
        SubstitutionVariant variant = (SubstitutionVariant) parser.parse(tatdn2, tatdn2Canonical, "c.*38+1G>A");
        assertEquals("G", variant.Ref);
        assertEquals("A", variant.Alt);
        assertEquals(38, ((InIntronDownstreamOfCodingEnd) variant.AddressOfChangeStart).IndexOfExonicBase);
        assertEquals(1, ((InIntronDownstreamOfCodingEnd) variant.AddressOfChangeStart).RelativePositionOfIntronicBase);

        variant = (SubstitutionVariant) parser.parse(tatdn2, tatdn2Canonical, "*39-3C>A");
        assertEquals("C", variant.Ref);
        assertEquals("A", variant.Alt);
        assertEquals(39, ((InIntronDownstreamOfCodingEnd) variant.AddressOfChangeStart).IndexOfExonicBase);
        assertEquals(-3, ((InIntronDownstreamOfCodingEnd) variant.AddressOfChangeStart).RelativePositionOfIntronicBase);
    }

    @Test
    public void parseSingleBaseDeletion()
    {
        DeletionVariant deletion = (DeletionVariant) parser.parse(tatdn2, tatdn2Canonical, "c.420delT");
        assertEquals(420, ((InExon)deletion.AddressOfChangeStart).IndexOfBaseInCodingBases);

        // The base is not needed (and HGVS says it's best left out)
        deletion = (DeletionVariant) parser.parse(tatdn2, tatdn2Canonical, "c.420del");
        assertEquals(420, ((InExon)deletion.AddressOfChangeStart).IndexOfBaseInCodingBases);

        // The c. is optional too
        deletion = (DeletionVariant) parser.parse(tatdn2, tatdn2Canonical, "420del");
        assertEquals(420, ((InExon)deletion.AddressOfChangeStart).IndexOfBaseInCodingBases);
    }

    @Test
    public void parseRangeDeletion()
    {
        DeletionVariant deletion = (DeletionVariant) parser.parse(zyx, zyxCanonical, "c.3_7delGCGGC");
        assertEquals(3, ((InExon)deletion.AddressOfChangeStart).IndexOfBaseInCodingBases);
        assertEquals(7, ((InExon)deletion.AddressOfChangeEnd).IndexOfBaseInCodingBases);
    }

    @Test
    public void parseRangeDuplication()
    {
        DuplicationVariant duplication = (DuplicationVariant) parser.parse(zyx, zyxCanonical, "c.3_7dupGCGGC");
        assertEquals(3, ((InExon)duplication.AddressOfChangeStart).IndexOfBaseInCodingBases);
        assertEquals(7, ((InExon)duplication.AddressOfChangeEnd).IndexOfBaseInCodingBases);
    }

    @Test
    public void parseInsertion()
    {
        InsertionVariant insertionVariant = (InsertionVariant) parser.parse(zyx, zyxCanonical, "c.3_7insGCGGC");
        assertEquals(3, ((InExon)insertionVariant.AddressOfChangeStart).IndexOfBaseInCodingBases);
        assertEquals(7, ((InExon)insertionVariant.AddressOfChangeEnd).IndexOfBaseInCodingBases);
        assertEquals("GCGGC", insertionVariant.insertedBases());
    }
}
