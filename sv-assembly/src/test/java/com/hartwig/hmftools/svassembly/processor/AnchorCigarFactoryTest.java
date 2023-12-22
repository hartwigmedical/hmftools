package com.hartwig.hmftools.svassembly.processor;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;

import com.hartwig.hmftools.svassembly.TestUtils;
import com.hartwig.hmftools.svassembly.models.Alignment;

import org.apache.commons.lang3.tuple.Pair;
import org.junit.Test;

import htsjdk.samtools.Cigar;

public class AnchorCigarFactoryTest
{
    private Pair<Pair<Cigar, Integer>, Pair<Cigar, Integer>> test(
            final List<Alignment> alignments, final int leftIndex, final int rightIndex)
    {
        int sequencePosition = 1;
        for(final Alignment alignment : alignments)
            if(sequencePosition != alignment.SequenceStartPosition)
                throw new IllegalStateException("Invalid Alignment list");
            else
                sequencePosition += alignment.Length;

        final AnchorCigarFactory factory = new AnchorCigarFactory(TestUtils.config());
        return factory.anchorCigar(alignments, alignments.get(leftIndex), alignments.get(rightIndex));
    }

    @Test
    public void fullMatchTranslocation()
    {
        final var result = test(List.of(
                new Alignment("1", 1, 1, 100, false, 60),
                new Alignment("2", 1, 101, 200, false, 60)
        ), 0, 1);

        final Pair<Cigar, Integer> left = result.getLeft();
        assertThat(left.getLeft().toString()).isEqualTo("100M");
        assertThat(left.getRight()).isEqualTo(100);

        final Pair<Cigar, Integer> right = result.getRight();
        assertThat(right.getLeft().toString()).isEqualTo("200M");
        assertThat(right.getRight()).isEqualTo(200);
    }

    @Test
    public void cutOffBySoftClip()
    {
        final var result = test(List.of(
                Alignment.unmapped(1, 50),
                new Alignment("1", 1, 51, 100, false, 60),
                Alignment.unmapped(151, 10),
                new Alignment("2", 1, 161, 200, false, 60),
                Alignment.unmapped(361, 50)
        ), 1, 3);

        final Pair<Cigar, Integer> left = result.getLeft();
        assertThat(left.getLeft().toString()).isEqualTo("50S100M");
        assertThat(left.getRight()).isEqualTo(100);

        final Pair<Cigar, Integer> right = result.getRight();
        assertThat(right.getLeft().toString()).isEqualTo("200M50S");
        assertThat(right.getRight()).isEqualTo(200);
    }

    @Test
    public void cutOffByLargeInsert()
    {
        final var result = test(List.of(
                new Alignment("1", 1, 1, 50, false, 60),
                Alignment.unmapped(51, 50),
                new Alignment("1", 51, 101, 50, false, 60),
                Alignment.unmapped(151, 10),
                new Alignment("2", 1, 161, 100, false, 60),
                Alignment.unmapped(261, 50),
                new Alignment("2", 101, 311, 100, false, 60)
        ), 2, 4);

        final Pair<Cigar, Integer> left = result.getLeft();
        assertThat(left.getLeft().toString()).isEqualTo("100S50M");
        assertThat(left.getRight()).isEqualTo(50);

        final Pair<Cigar, Integer> right = result.getRight();
        assertThat(right.getLeft().toString()).isEqualTo("100M150S");
        assertThat(right.getRight()).isEqualTo(100);
    }

    @Test
    public void cutOffByLargeInsertInvert()
    {
        final var result = test(List.of(
                new Alignment("1", 51, 1, 50, true, 60),
                Alignment.unmapped(51, 50),
                new Alignment("1", 1, 101, 50, true, 60),
                Alignment.unmapped(151, 10),
                new Alignment("2", 101, 161, 100, true, 60),
                Alignment.unmapped(261, 50),
                new Alignment("2", 1, 311, 100, true, 60)
        ), 2, 4);

        final Pair<Cigar, Integer> left = result.getLeft();
        assertThat(left.getLeft().toString()).isEqualTo("50M100S");
        assertThat(left.getRight()).isEqualTo(50);

        final Pair<Cigar, Integer> right = result.getRight();
        assertThat(right.getLeft().toString()).isEqualTo("150S100M");
        assertThat(right.getRight()).isEqualTo(100);
    }

    @Test
    public void cutOffByLargeDelete()
    {
        final var result = test(List.of(
                new Alignment("1", 1, 1, 50, false, 60),
                new Alignment("1", 151, 51, 50, false, 60),
                Alignment.unmapped(101, 10),
                new Alignment("2", 1, 111, 100, false, 60),
                new Alignment("2", 201, 211, 100, false, 60)
        ), 1, 3);

        final Pair<Cigar, Integer> left = result.getLeft();
        assertThat(left.getLeft().toString()).isEqualTo("50S50M");
        assertThat(left.getRight()).isEqualTo(50);

        final Pair<Cigar, Integer> right = result.getRight();
        assertThat(right.getLeft().toString()).isEqualTo("100M100S");
        assertThat(right.getRight()).isEqualTo(100);
    }

    @Test
    public void notCutOffBySmallInserts()
    {
        final var result = test(List.of(
                new Alignment("1", 1, 1, 50, false, 60),
                Alignment.unmapped(51, 10),
                new Alignment("1", 51, 61, 50, false, 60),
                Alignment.unmapped(111, 10),
                new Alignment("2", 1, 121, 100, false, 60),
                Alignment.unmapped(221, 10),
                new Alignment("2", 101, 231, 100, false, 60)
        ), 2, 4);

        final Pair<Cigar, Integer> left = result.getLeft();
        assertThat(left.getLeft().toString()).isEqualTo("50M10I50M");
        assertThat(left.getRight()).isEqualTo(100);

        final Pair<Cigar, Integer> right = result.getRight();
        assertThat(right.getLeft().toString()).isEqualTo("100M10I100M");
        assertThat(right.getRight()).isEqualTo(200);
    }

    @Test
    public void notCutOffBySmallDeletes()
    {
        final var result = test(List.of(
                new Alignment("1", 1, 1, 50, false, 60),
                new Alignment("1", 61, 51, 50, false, 60),
                Alignment.unmapped(101, 10),
                new Alignment("2", 1, 111, 100, false, 60),
                new Alignment("2", 111, 211, 100, false, 60)
        ), 1, 3);

        final Pair<Cigar, Integer> left = result.getLeft();
        assertThat(left.getLeft().toString()).isEqualTo("50M10D50M");
        assertThat(left.getRight()).isEqualTo(110);

        final Pair<Cigar, Integer> right = result.getRight();
        assertThat(right.getLeft().toString()).isEqualTo("100M10D100M");
        assertThat(right.getRight()).isEqualTo(210);
    }
}