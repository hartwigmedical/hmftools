package com.hartwig.hmftools.orange.algo.protect;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectTestFactory;

import org.junit.Test;

public class EvidenceSelectorTest {

    @Test
    public void canSelectEventsWithEvidence() {
        List<ProtectEvidence> evidences = Lists.newArrayList();
        evidences.add(ProtectTestFactory.builder().gene("gene A").event("event A").build());
        evidences.add(ProtectTestFactory.builder().gene(null).event("event B").build());

        assertTrue(EvidenceSelector.hasEvidence(evidences, "gene A", "event A"));
        assertTrue(EvidenceSelector.hasEvidence(evidences, null, "event B"));

        assertFalse(EvidenceSelector.hasEvidence(evidences, "gene A", "event B"));
        assertFalse(EvidenceSelector.hasEvidence(evidences, "gene B", "event B"));
    }
}