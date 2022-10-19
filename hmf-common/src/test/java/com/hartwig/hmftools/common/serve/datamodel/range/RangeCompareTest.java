package com.hartwig.hmftools.common.serve.datamodel.range;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.serve.datamodel.MutationTypeFilter;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class RangeCompareTest {

    @Test
    public void canCompareRanges() {
        RangeAnnotation annotation1 = create("X", 1, 10, "gene", "transcript", 1, MutationTypeFilter.INFRAME_DELETION);
        RangeAnnotation annotation2 = create("X", 1, 10, "gene", "transcript", 1, MutationTypeFilter.INFRAME_INSERTION);
        RangeAnnotation annotation3 = create("X", 1, 10, "a other gene", "transcript", 1, MutationTypeFilter.INFRAME_DELETION);
        RangeAnnotation annotation4 = create("X", 2, 10, "gene", "transcript", 1, MutationTypeFilter.INFRAME_DELETION);

        assertTrue(RangeCompare.compare(annotation1, annotation2) < 0);
        assertTrue(RangeCompare.compare(annotation1, annotation3) > 0);
        assertTrue(RangeCompare.compare(annotation4, annotation1) > 0);

        assertEquals(0, RangeCompare.compare(annotation1, annotation1));
    }

    @NotNull
    private static RangeAnnotation create(@NotNull String chromosome, int start, int end, @NotNull String gene, @NotNull String transcript,
            int rank, @NotNull MutationTypeFilter mutationTypeFilter) {
        return new RangeAnnotation() {
            @NotNull
            @Override
            public String gene() {
                return gene;
            }

            @NotNull
            @Override
            public String transcript() {
                return transcript;
            }

            @Override
            public int rank() {
                return rank;
            }

            @NotNull
            @Override
            public MutationTypeFilter mutationType() {
                return mutationTypeFilter;
            }

            @NotNull
            @Override
            public String chromosome() {
                return chromosome;
            }

            @Override
            public int start() {
                return start;
            }

            @Override
            public int end() {
                return end;
            }
        };
    }
}