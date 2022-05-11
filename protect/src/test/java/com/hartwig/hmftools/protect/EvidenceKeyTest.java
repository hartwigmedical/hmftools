package com.hartwig.hmftools.protect;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectTestFactory;

import org.junit.Test;

public class EvidenceKeyTest {

    @Test
    public void canCreateUniqueEventSet() {
        List<ProtectEvidence> evidences = Lists.newArrayList();

        ProtectEvidence evidence1 = ProtectTestFactory.builder().gene("gene 1").event("event 1").treatment("treatment 1").build();
        evidences.add(evidence1);
        evidences.add(evidence1);
        evidences.add(evidence1);

        Set<EvidenceKey> keys = EvidenceKey.buildUniqueEventSet(evidences);
        assertEquals(1, keys.size());
        EvidenceKey key = keys.iterator().next();
        assertEquals("gene 1", key.gene());
        assertEquals("event 1", key.event());
        assertEquals("treatment 1", key.treatment());

        ProtectEvidence evidence2 = ProtectTestFactory.builder().gene(null).event("event 2").treatment("treatment 2").build();
        evidences.add(evidence2);
        evidences.add(evidence2);
        evidences.add(evidence2);

        assertEquals(2, EvidenceKey.buildUniqueEventSet(evidences).size());
    }
}