package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.redux.consensus.NonStandardBaseBuilder.alignAnnotatedReads;

import static org.junit.Assert.assertEquals;

import static htsjdk.samtools.CigarOperator.M;

import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.redux.consensus.NonStandardBaseBuilder.AnnotatedBase;
import com.hartwig.hmftools.redux.consensus.NonStandardBaseBuilder.ExtendedRefPos;

import org.junit.Test;

public class NonStandardBaseBuilderTest
{
    @Test
    public void testAlignAnnotatedReadsDelSkip()
    {
        List<AnnotatedBase> annotatedRead1 = Lists.newArrayList(
                new AnnotatedBase(new ExtendedRefPos(1, 0), (byte) 'A', (byte) 1, M),
                new AnnotatedBase(new ExtendedRefPos(3, 0), (byte) 'A', (byte) 1, M));

        List<AnnotatedBase> annotatedRead2 = Lists.newArrayList(
                new AnnotatedBase(new ExtendedRefPos(1, 0), (byte) 'A', (byte) 1, M),
                new AnnotatedBase(new ExtendedRefPos(2, 0), (byte) 'A', (byte) 1, M),
                new AnnotatedBase(new ExtendedRefPos(3, 0), (byte) 'A', (byte) 1, M)
        );

        List<List<AnnotatedBase>> annotatedReads = List.of(annotatedRead1, annotatedRead2);
        Collection<List<AnnotatedBase>> alignment = alignAnnotatedReads(annotatedReads);
        List<Integer> actualDepth = alignment.stream().map(List::size).collect(Collectors.toList());
        List<Integer> expectedDepth = List.of(2, 1, 2);

        assertEquals(expectedDepth, actualDepth);
    }
}
