package com.hartwig.hmftools.sage.testutil;

import static com.hartwig.hmftools.sage.testutil.MutatedBases.getAlignedBases;
import static com.hartwig.hmftools.sage.testutil.MutatedBases.getCigarStr;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import static htsjdk.samtools.CigarOperator.EQ;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.X;

import java.util.List;

import com.google.common.collect.Lists;

import org.junit.Test;

public class MutatesBasesTest
{
    @Test
    public void testGetCigarStrEmpty()
    {
        assertEquals("", getCigarStr(Lists.newArrayList()));
    }

    @Test
    public void testGetCigarStrNoDel()
    {
        List<MutatedBases.MutatedBase> bases = Lists.newArrayList(
                new MutatedBases.MutatedBase(1, 'A', EQ),
                new MutatedBases.MutatedBase(2, 'A', X),
                new MutatedBases.MutatedBase(2, 'A', I),
                new MutatedBases.MutatedBase(3, 'A', EQ)
        );
        String expected = "2M1I1M";
        String actual = getCigarStr(bases);

        assertEquals(expected, actual);
    }

    @Test
    public void testGetCigarStrWithDel()
    {
        List<MutatedBases.MutatedBase> bases = Lists.newArrayList(
                new MutatedBases.MutatedBase(1, 'A', EQ),
                new MutatedBases.MutatedBase(2, 'A', X),
                new MutatedBases.MutatedBase(2, 'A', I),
                new MutatedBases.MutatedBase(6, 'A', EQ)
        );
        String expected = "2M1I3D1M";
        String actual = getCigarStr(bases);

        assertEquals(expected, actual);
    }

    @Test
    public void testGetCigarStrAlignedBasesLeftSoftClip()
    {
        List<MutatedBases.MutatedBase> bases = Lists.newArrayList(
                new MutatedBases.MutatedBase(1, 'A', EQ),
                new MutatedBases.MutatedBase(2, 'A', X),
                new MutatedBases.MutatedBase(2, 'A', I),
                new MutatedBases.MutatedBase(6, 'A', EQ)
        );
        MutatedBases.AlignedBases aligned = new MutatedBases.AlignedBases(5, 0, bases);
        String expected = "5S2M1I3D1M";
        String actual = getCigarStr(aligned);

        assertEquals(expected, actual);
    }

    @Test
    public void testGetAlignedBasesEmpty()
    {
        assertNull(getAlignedBases(null, Lists.newArrayList()));
    }

    @Test
    public void testGetAlignedBasesAllAligned()
    {
        String refBases = null;
        List<MutatedBases.MutatedBase> bases = Lists.newArrayList(
                new MutatedBases.MutatedBase(1, 'A', EQ),
                new MutatedBases.MutatedBase(2, 'A', X),
                new MutatedBases.MutatedBase(2, 'A', I),
                new MutatedBases.MutatedBase(6, 'A', EQ)
        );
        MutatedBases.AlignedBases expected = new MutatedBases.AlignedBases(0, 0, bases);
        MutatedBases.AlignedBases actual = getAlignedBases(refBases, bases);

        assertEquals(expected, actual);
    }

    @Test
    public void testGetAlignedBasesAbsorbedEndsNoSoftClip()
    {
        String refBases = "GGGAAAAAAA";
        List<MutatedBases.MutatedBase> bases = Lists.newArrayList(
                new MutatedBases.MutatedBase(1, 'A', X),
                new MutatedBases.MutatedBase(1, 'A', I),
                new MutatedBases.MutatedBase(2, 'A', X),
                new MutatedBases.MutatedBase(7, 'A', EQ),
                new MutatedBases.MutatedBase(8, 'A', X),
                new MutatedBases.MutatedBase(8, 'A', I),
                new MutatedBases.MutatedBase(9, 'A', EQ),
                new MutatedBases.MutatedBase(9, 'A', I)
        );
        List<MutatedBases.MutatedBase> expectedBases = Lists.newArrayList(
                new MutatedBases.MutatedBase(4, 'A', EQ),
                new MutatedBases.MutatedBase(5, 'A', EQ),
                new MutatedBases.MutatedBase(6, 'A', EQ),
                new MutatedBases.MutatedBase(7, 'A', EQ),
                new MutatedBases.MutatedBase(8, 'A', X),
                new MutatedBases.MutatedBase(8, 'A', I),
                new MutatedBases.MutatedBase(9, 'A', EQ),
                new MutatedBases.MutatedBase(10, 'A', EQ)
        );
        MutatedBases.AlignedBases expected = new MutatedBases.AlignedBases(0, 0, expectedBases);
        MutatedBases.AlignedBases actual = getAlignedBases(refBases, bases);

        assertEquals(expected, actual);
    }

    @Test
    public void testGetAlignedBasesSoftClip()
    {
        String refBases = "GGGGGGAAAG";
        List<MutatedBases.MutatedBase> bases = Lists.newArrayList(
                new MutatedBases.MutatedBase(1, 'A', X),
                new MutatedBases.MutatedBase(1, 'A', I),
                new MutatedBases.MutatedBase(2, 'A', X),
                new MutatedBases.MutatedBase(7, 'A', EQ),
                new MutatedBases.MutatedBase(8, 'A', X),
                new MutatedBases.MutatedBase(8, 'A', I),
                new MutatedBases.MutatedBase(9, 'A', EQ),
                new MutatedBases.MutatedBase(9, 'A', I)
        );
        List<MutatedBases.MutatedBase> expectedBases = Lists.newArrayList(
                new MutatedBases.MutatedBase(7, 'A', EQ),
                new MutatedBases.MutatedBase(8, 'A', X),
                new MutatedBases.MutatedBase(8, 'A', I),
                new MutatedBases.MutatedBase(9, 'A', EQ)
        );
        MutatedBases.AlignedBases expected = new MutatedBases.AlignedBases(3, 1, expectedBases);
        MutatedBases.AlignedBases actual = getAlignedBases(refBases, bases);

        assertEquals(expected, actual);
    }

    @Test
    public void testGetAlignedBasesEmptyAlignment()
    {
        String refBases = null;
        List<MutatedBases.MutatedBase> bases = Lists.newArrayList(
                new MutatedBases.MutatedBase(1, 'A', X),
                new MutatedBases.MutatedBase(2, 'A', X),
                new MutatedBases.MutatedBase(3, 'A', X)
        );

        assertNull(getAlignedBases(refBases, bases));
    }

    @Test
    public void testLeftMutIndexFromRefPosNoBases()
    {
        MutatedBases bases = new MutatedBases(null, Lists.newArrayList());
        assertEquals(-1, bases.leftMutIndexFromRefPos(2));
    }

    @Test
    public void testLeftMutIndexFromRefPosNoIndex()
    {
        String refBases = "AAAA";
        MutatedBasesBuilder basesBuilder = new MutatedBasesBuilder(refBases);
        basesBuilder.delBases(1, 2);

        MutatedBases bases = basesBuilder.build();
        assertEquals(-1, bases.leftMutIndexFromRefPos(2));
    }

    @Test
    public void testLeftMutIndexFromRefPosExactMatch()
    {
        String refBases = "AAAA";
        MutatedBasesBuilder basesBuilder = new MutatedBasesBuilder(refBases);
        MutatedBases bases = basesBuilder.build();
        assertEquals(1, bases.leftMutIndexFromRefPos(2));
    }

    @Test
    public void testLeftMutIndexFromRefPosNonExactMatch()
    {
        String refBases = "AAAA";
        MutatedBasesBuilder basesBuilder = new MutatedBasesBuilder(refBases);
        basesBuilder.delBases(2, 1);

        MutatedBases bases = basesBuilder.build();
        assertEquals(0, bases.leftMutIndexFromRefPos(2));
    }

    @Test
    public void testRightMutIndexFromRefPosNoBases()
    {
        MutatedBases bases = new MutatedBases(null, Lists.newArrayList());
        assertEquals(-1, bases.rightMutIndexFromRefPos(2));
    }

    @Test
    public void testRightMutIndexFromRefPosNoIndex()
    {
        String refBases = "AAAA";
        MutatedBasesBuilder basesBuilder = new MutatedBasesBuilder(refBases);
        basesBuilder.delBases(3, 2);

        MutatedBases bases = basesBuilder.build();
        assertEquals(-1, bases.rightMutIndexFromRefPos(3));
    }

    @Test
    public void testRightMutIndexFromRefPosExactMatch()
    {
        String refBases = "AAAA";
        MutatedBasesBuilder basesBuilder = new MutatedBasesBuilder(refBases);
        MutatedBases bases = basesBuilder.build();
        assertEquals(2, bases.rightMutIndexFromRefPos(3));
    }

    @Test
    public void testRightMutIndexFromRefPosNonExactMatch()
    {
        String refBases = "AAAA";
        MutatedBasesBuilder basesBuilder = new MutatedBasesBuilder(refBases);
        basesBuilder.delBases(3, 1);

        MutatedBases bases = basesBuilder.build();
        assertEquals(2, bases.rightMutIndexFromRefPos(3));
    }
}
