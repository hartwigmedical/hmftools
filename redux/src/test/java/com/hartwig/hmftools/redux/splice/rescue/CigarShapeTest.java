package com.hartwig.hmftools.redux.splice.rescue;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;

public class CigarShapeTest
{
    @Test
    public void testParseSimple()
    {
        final List<CigarShape.Element> elements = CigarShape.parse("94M57S");
        assertEquals(2, elements.size());
        assertEquals(94, elements.get(0).Length);
        assertEquals('M', elements.get(0).Op);
        assertEquals(57, elements.get(1).Length);
        assertEquals('S', elements.get(1).Op);
    }

    @Test
    public void testParseMultipleOps()
    {
        final List<CigarShape.Element> elements = CigarShape.parse("50M200N40M61S");
        assertEquals(4, elements.size());
        assertEquals('N', elements.get(1).Op);
        assertEquals(200, elements.get(1).Length);
    }

    @Test
    public void testParseEmpty()
    {
        assertTrue(CigarShape.parse("").isEmpty());
        assertTrue(CigarShape.parse(null).isEmpty());
        assertTrue(CigarShape.parse("*").isEmpty());
    }

    @Test(expected = IllegalArgumentException.class)
    public void testParseInvalidCigar()
    {
        CigarShape.parse("94MXX57S");
    }

    @Test
    public void testFormat()
    {
        final String cigar = "94M79N57M";
        assertEquals(cigar, CigarShape.format(CigarShape.parse(cigar)));
    }

    @Test
    public void testReadLength()
    {
        assertEquals(151, CigarShape.readLength(CigarShape.parse("94M57S")));
        assertEquals(151, CigarShape.readLength(CigarShape.parse("94M79N57M")));
        assertEquals(154, CigarShape.readLength(CigarShape.parse("94M3I57S")));
        assertEquals(151, CigarShape.readLength(CigarShape.parse("94M3D57S")));
    }

    @Test
    public void testReferenceSpan()
    {
        assertEquals(94, CigarShape.referenceSpan(CigarShape.parse("94M57S")));
        assertEquals(94 + 79 + 57, CigarShape.referenceSpan(CigarShape.parse("94M79N57M")));
        assertEquals(94, CigarShape.referenceSpan(CigarShape.parse("94M3I57S")));
        assertEquals(94 + 3, CigarShape.referenceSpan(CigarShape.parse("94M3D54S")));
    }

    @Test
    public void testLeadingTrailingSoftClip()
    {
        assertEquals(0, CigarShape.leadingSoftClip(CigarShape.parse("94M57S")));
        assertEquals(57, CigarShape.trailingSoftClip(CigarShape.parse("94M57S")));
        assertEquals(57, CigarShape.leadingSoftClip(CigarShape.parse("57S94M")));
        assertEquals(0, CigarShape.trailingSoftClip(CigarShape.parse("57S94M")));
        assertEquals(10, CigarShape.leadingSoftClip(CigarShape.parse("10S94M47S")));
        assertEquals(47, CigarShape.trailingSoftClip(CigarShape.parse("10S94M47S")));
    }

    @Test
    public void testHasHardClip()
    {
        assertFalse(CigarShape.hasHardClip(CigarShape.parse("94M57S")));
        assertTrue(CigarShape.hasHardClip(CigarShape.parse("10H94M57S")));
        assertTrue(CigarShape.hasHardClip(CigarShape.parse("94M57S10H")));
    }

    @Test
    public void testMatchedBases()
    {
        assertEquals(94, CigarShape.matchedBases(CigarShape.parse("94M57S")));
        assertEquals(94 + 57, CigarShape.matchedBases(CigarShape.parse("94M79N57M")));
        assertEquals(94 + 54, CigarShape.matchedBases(CigarShape.parse("94M3I54M")));
    }
}
