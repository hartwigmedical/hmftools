package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.Utils.getBestScoringElement;
import static com.hartwig.hmftools.panelbuilder.Utils.outwardMovingOffsets;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Optional;
import java.util.stream.Stream;

import org.junit.Test;

public class UtilsTest
{
    @Test
    public void testGetBestScoringElementEmptyStream()
    {
        assertEquals(Optional.empty(), Utils.<Double>getBestScoringElement(Stream.empty(), e -> e, s -> s == 0, true));

        assertEquals(Optional.empty(), Utils.<Double>getBestScoringElement(Stream.empty(), e -> e, s -> s == 0, false));
    }

    @Test
    public void testGetBestScoringElementOptimal()
    {
        assertEquals(Optional.of(3), getBestScoringElement(Stream.of(1, 2, 3, 4, 5), e -> e, s -> s == 3, true));
        assertEquals(Optional.of(3), getBestScoringElement(Stream.of(3, 2, 1, 4, 5), e -> e, s -> s == 3, true));

        assertEquals(Optional.of(3), getBestScoringElement(Stream.of(5, 4, 3, 2, 1), e -> e, s -> s == 3, false));
        assertEquals(Optional.of(3), getBestScoringElement(Stream.of(3, 4, 5, 2, 1), e -> e, s -> s == 3, false));
    }

    @Test
    public void testGetBestScoringElementNoOptimal()
    {
        assertEquals(Optional.of(6), getBestScoringElement(Stream.of(1, 2, 6, 4, 5), e -> e, s -> s == 0, true));
        assertEquals(Optional.of(6), getBestScoringElement(Stream.of(6, 2, 1, 4, 5), e -> e, s -> s == 0, true));

        assertEquals(Optional.of(0), getBestScoringElement(Stream.of(5, 0, 3, 2, 1), e -> e, s -> s == 0, false));
        assertEquals(Optional.of(0), getBestScoringElement(Stream.of(0, 5, 3, 2, 1), e -> e, s -> s == 0, false));
    }

    @Test
    public void testGetBestScoringElementScoreFunc()
    {
        assertEquals(Optional.of(6), getBestScoringElement(Stream.of(2, 4, 6, 8, 10), e -> (double) e / 2, s -> s == 3.0, true));
    }

    @Test
    public void testOutwardMovingOffsets()
    {
        assertEquals(List.of(0), outwardMovingOffsets(0, 0).boxed().toList());
        assertEquals(List.of(0, 1, -1, 2, -2, 3, -3), outwardMovingOffsets(-3, 3).boxed().toList());
        assertEquals(List.of(0, 1, -1, 2, -2, 3, -3, 4, -4), outwardMovingOffsets(-4, 4).boxed().toList());
        assertEquals(List.of(0, 1, 2, 3), outwardMovingOffsets(0, 3).boxed().toList());
        assertEquals(List.of(0, -1, -2, -3), outwardMovingOffsets(-3, 0).boxed().toList());
        assertEquals(List.of(0, 1, -1, 2, -2, 3, -3, 4, -4, 5, 6, 7), outwardMovingOffsets(-4, 7).boxed().toList());
    }
}
