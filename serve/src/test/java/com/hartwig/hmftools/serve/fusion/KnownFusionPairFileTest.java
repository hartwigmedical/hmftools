package com.hartwig.hmftools.serve.fusion;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.apache.commons.compress.utils.Lists;
import org.junit.Test;

public class KnownFusionPairFileTest {

    @Test
    public void canSortFusionPairs() {
        KnownFusionPair pair1  = ImmutableKnownFusionPair.builder().geneUp("A").geneDown("B").minExonUp(2).minExonDown(3).build();
        KnownFusionPair pair2 = ImmutableKnownFusionPair.builder().geneUp("A").geneDown("B").minExonUp(5).minExonDown(6).build();
        KnownFusionPair pair3 = ImmutableKnownFusionPair.builder().geneUp("A").geneDown("B").build();
        KnownFusionPair pair4 = ImmutableKnownFusionPair.builder().geneUp("A").geneDown("C").build();
        KnownFusionPair pair5 = ImmutableKnownFusionPair.builder().geneUp("X").geneDown("A").build();

        List<KnownFusionPair> knownFusionPairs = Lists.newArrayList();
        knownFusionPairs.add(pair3);
        knownFusionPairs.add(pair1);
        knownFusionPairs.add(pair5);
        knownFusionPairs.add(pair4);
        knownFusionPairs.add(pair2);

        List<KnownFusionPair> sorted = KnownFusionPairFile.sort(knownFusionPairs);

        assertEquals(pair1, sorted.get(0));
        assertEquals(pair2, sorted.get(1));
        assertEquals(pair3, sorted.get(2));
        assertEquals(pair4, sorted.get(3));
        assertEquals(pair5, sorted.get(4));
    }

}