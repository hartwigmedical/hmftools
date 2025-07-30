package com.hartwig.hmftools.geneutils.probequality;

import static com.hartwig.hmftools.geneutils.probequality.Utils.partitionStream;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.junit.Test;

public class UtilsTest
{
    @Test
    public void testPartitionStreamSingle()
    {
        List<Integer> regions = List.of(1, 2);
        List<List<Integer>> expected = List.of(List.of(1, 2));
        List<List<Integer>> actual = partitionStream(regions.stream(), 3).toList();
        assertEquals(expected, actual);
    }

    @Test
    public void testPartitionStreamMultiple()
    {
        List<Integer> regions = List.of(1, 2, 3, 4, 5, 6, 7);
        List<List<Integer>> expected = List.of(
                List.of(1, 2, 3),
                List.of(4, 5, 6),
                List.of(7)
        );
        List<List<Integer>> actual = partitionStream(regions.stream(), 3).toList();
        assertEquals(expected, actual);
    }
}
