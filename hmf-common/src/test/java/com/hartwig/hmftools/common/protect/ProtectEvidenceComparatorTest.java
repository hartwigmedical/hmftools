package com.hartwig.hmftools.common.protect;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;

import org.junit.Test;

public class ProtectEvidenceComparatorTest {

    @Test
    public void canSortProtectEvidence() {
        ProtectEvidence evidence1 = ProtectTestFactory.builder().reported(true).gene("A").build();
        ProtectEvidence evidence2 = ProtectTestFactory.builder().reported(false).gene("A").build();
        ProtectEvidence evidence3 = ProtectTestFactory.builder().reported(true).gene("B").build();

        List<ProtectEvidence> evidences = Lists.newArrayList(evidence1, evidence2, evidence3);
        evidences.sort(new EvidenceComparator());

        assertEquals(evidence1, evidences.get(0));
        assertEquals(evidence3, evidences.get(1));
        assertEquals(evidence2, evidences.get(2));
    }
}