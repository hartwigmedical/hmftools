package com.hartwig.hmftools.sage.testutil;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import static htsjdk.samtools.CigarOperator.EQ;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.X;

import java.util.List;

import com.beust.jcommander.internal.Lists;
import com.hartwig.hmftools.common.region.BaseRegion;

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

    @Test
    public void testMutationBounds()
    {
        String refBases = "ATGC".repeat(2);

        // no mutation
        MutatedBasesBuilder builder = new MutatedBasesBuilder(refBases);
        assertNull(builder.build().refPosMutationBounds());

        // mutate base to ref
        builder = new MutatedBasesBuilder(refBases);
        builder.mutateBase(3, 'G');
        assertNull(builder.build().refPosMutationBounds());

        // mutate base back to ref
        builder = new MutatedBasesBuilder(refBases);
        builder.mutateBase(3, 'T');
        builder.mutateBase(3, 'G');
        assertNull(builder.build().refPosMutationBounds());

        // snv
        builder = new MutatedBasesBuilder(refBases);
        builder.mutateBase(3, 'T');
        assertEquals(new BaseRegion(3, 3), builder.build().refPosMutationBounds());

        // insert
        builder = new MutatedBasesBuilder(refBases);
        builder.insertBases(3, "AAA");
        assertEquals(new BaseRegion(3, 4), builder.build().refPosMutationBounds());

        // del
        builder = new MutatedBasesBuilder(refBases);
        builder.delBases(3, 4);
        assertEquals(new BaseRegion(3, 6), builder.build().refPosMutationBounds());

        // combined
        builder = new MutatedBasesBuilder(refBases);
        builder.mutateBase(8, 'C');
        builder.mutateBase(4, 'T');
        builder.insertBases(1, "AA");
        builder.delBases(5, 2);
        assertEquals(new BaseRegion(1, 6), builder.build().refPosMutationBounds());
    }
}
