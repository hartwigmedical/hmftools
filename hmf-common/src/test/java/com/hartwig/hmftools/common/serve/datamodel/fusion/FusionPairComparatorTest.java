package com.hartwig.hmftools.common.serve.datamodel.fusion;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class FusionPairComparatorTest {

    @Test
    public void canSortFusionPairs() {
        FusionPair pair1 = create("A", "B", 2, 3);
        FusionPair pair2 = create("A", "B", null, null);
        FusionPair pair3 = create("A", "C", null, null);
        FusionPair pair4 = create("B", "C", null, null);
        FusionPair pair5 = create("X", "A", null, null);

        List<FusionPair> fusionPairs = Lists.newArrayList(pair4, pair1, pair5, pair3, pair2);
        fusionPairs.sort(new FusionPairComparator());

        assertEquals(pair1, fusionPairs.get(0));
        assertEquals(pair2, fusionPairs.get(1));
        assertEquals(pair3, fusionPairs.get(2));
        assertEquals(pair4, fusionPairs.get(3));
        assertEquals(pair5, fusionPairs.get(4));
    }

    @NotNull
    private static FusionPair create(@NotNull String geneUp, @NotNull String geneDown, @Nullable Integer minExonUp,
            @Nullable Integer minExonDown) {
        return new FusionPair() {
            @NotNull
            @Override
            public String geneUp() {
                return geneUp;
            }

            @Nullable
            @Override
            public Integer minExonUp() {
                return minExonUp;
            }

            @Nullable
            @Override
            public Integer maxExonUp() {
                return null;
            }

            @NotNull
            @Override
            public String geneDown() {
                return geneDown;
            }

            @Nullable
            @Override
            public Integer minExonDown() {
                return minExonDown;
            }

            @Nullable
            @Override
            public Integer maxExonDown() {
                return null;
            }
        };
    }
}