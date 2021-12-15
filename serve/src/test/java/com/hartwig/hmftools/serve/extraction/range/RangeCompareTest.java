package com.hartwig.hmftools.serve.extraction.range;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.serve.extraction.exon.ImmutableExonAnnotation;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilter;

import org.junit.Test;

public class RangeCompareTest {

    @Test
    public void canCompareRanges() {
        RangeAnnotation annotation1 = ImmutableExonAnnotation.builder()
                .chromosome("X")
                .start(1)
                .end(10)
                .gene("gene")
                .transcript("transcript")
                .rank(1)
                .mutationType(MutationTypeFilter.INFRAME_DELETION)
                .build();

        RangeAnnotation annotation2 =
                ImmutableExonAnnotation.builder().from(annotation1).mutationType(MutationTypeFilter.INFRAME_INSERTION).build();

        RangeAnnotation annotation3 =
                ImmutableExonAnnotation.builder().from(annotation1).gene("a other gene").build();

        RangeAnnotation annotation4 =
                ImmutableExonAnnotation.builder().from(annotation1).start(2).build();

        assertTrue(RangeCompare.compare(annotation1, annotation2) < 0);
        assertTrue(RangeCompare.compare(annotation1, annotation3) > 0);
        assertTrue(RangeCompare.compare(annotation4, annotation1) > 0);

        assertEquals(0, RangeCompare.compare(annotation1, annotation1));
    }
}