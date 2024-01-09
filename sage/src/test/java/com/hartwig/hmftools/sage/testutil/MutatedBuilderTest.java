package com.hartwig.hmftools.sage.testutil;

import static org.junit.Assert.assertEquals;

import static htsjdk.samtools.CigarOperator.EQ;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.X;

import java.util.List;

import com.beust.jcommander.internal.Lists;

import org.junit.Test;

public class MutatedBuilderTest
{
    @Test
    public void testNoMutations()
    {
        String refBases = "ATGC";
        MutatedBasesBuilder builder = new MutatedBasesBuilder(refBases);
        List<MutatedBases.MutatedBase> actual = builder.build().mutatedBases();

        List<MutatedBases.MutatedBase> expected = Lists.newArrayList();
        for(int i = 0; i < refBases.length(); ++i)
        {
            expected.add(new MutatedBases.MutatedBase(i + 1, refBases.charAt(i), EQ));
        }

        assertEquals(expected, actual);
    }

    @Test
    public void testSNV()
    {
        String refBases = "ATGC";
        MutatedBasesBuilder builder = new MutatedBasesBuilder(refBases);
        builder.mutateBase(2, 'A');
        List<MutatedBases.MutatedBase> actual = builder.build().mutatedBases();
        List<MutatedBases.MutatedBase> expected = Lists.newArrayList(
                new MutatedBases.MutatedBase(1, 'A', EQ),
                new MutatedBases.MutatedBase(2, 'A', X),
                new MutatedBases.MutatedBase(3, 'G', EQ),
                new MutatedBases.MutatedBase(4, 'C', EQ)
        );

        assertEquals(expected, actual);
    }

    @Test
    public void testInsert()
    {
        String refBases = "ATGC";
        MutatedBasesBuilder builder = new MutatedBasesBuilder(refBases);
        builder.insertBases(2, "AA");
        List<MutatedBases.MutatedBase> actual = builder.build().mutatedBases();
        List<MutatedBases.MutatedBase> expected = Lists.newArrayList(
                new MutatedBases.MutatedBase(1, 'A', EQ),
                new MutatedBases.MutatedBase(2, 'T', EQ),
                new MutatedBases.MutatedBase(2, 'A', I),
                new MutatedBases.MutatedBase(2, 'A', I),
                new MutatedBases.MutatedBase(3, 'G', EQ),
                new MutatedBases.MutatedBase(4, 'C', EQ)
        );

        assertEquals(expected, actual);
    }

    @Test
    public void testDel()
    {
        String refBases = "ATGC";
        MutatedBasesBuilder builder = new MutatedBasesBuilder(refBases);
        builder.delBases(2, 2);
        List<MutatedBases.MutatedBase> actual = builder.build().mutatedBases();
        List<MutatedBases.MutatedBase> expected = Lists.newArrayList(
                new MutatedBases.MutatedBase(1, 'A', EQ),
                new MutatedBases.MutatedBase(4, 'C', EQ)
        );

        assertEquals(expected, actual);
    }
}
