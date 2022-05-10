package com.hartwig.hmftools.serve.extraction.fusion;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;

import org.junit.Test;

public class FusionPairComparatorTest {

    @Test
    public void canSortFusionPairs() {
        KnownFusionPair pair1 = ImmutableKnownFusionPair.builder().geneUp("A").geneDown("B").minExonUp(2).minExonDown(3).build();
        KnownFusionPair pair2 = ImmutableKnownFusionPair.builder().geneUp("A").geneDown("B").minExonUp(5).minExonDown(6).build();
        KnownFusionPair pair3 = ImmutableKnownFusionPair.builder().geneUp("A").geneDown("B").build();
        KnownFusionPair pair4 = ImmutableKnownFusionPair.builder().geneUp("A").geneDown("C").build();
        KnownFusionPair pair5 = ImmutableKnownFusionPair.builder().geneUp("X").geneDown("A").build();

        List<KnownFusionPair> knownFusionPairs = Lists.newArrayList(pair3, pair1, pair5, pair4, pair2);
        knownFusionPairs.sort(new FusionPairComparator());

        assertEquals(pair1, knownFusionPairs.get(0));
        assertEquals(pair2, knownFusionPairs.get(1));
        assertEquals(pair3, knownFusionPairs.get(2));
        assertEquals(pair4, knownFusionPairs.get(3));
        assertEquals(pair5, knownFusionPairs.get(4));
    }
}