package com.hartwig.hmftools.protect;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectTestFactory;

import org.junit.Test;

public class EventKeyTest {

    @Test
    public void canCreateUniqueEventSet() {
        List<ProtectEvidence> evidences = Lists.newArrayList();

        ProtectEvidence evidence1 =
                ProtectTestFactory.builder().gene("Gene1").transcript("123").isCanonical(true).event("Event 1").build();
        evidences.add(evidence1);
        evidences.add(evidence1);
        evidences.add(evidence1);

        Set<EventKey> keys = EventKey.buildUniqueEventSet(evidences);
        assertEquals(1, keys.size());
        EventKey key = keys.iterator().next();
        assertEquals("Gene1", key.gene());
        assertEquals("Event 1", key.event());

        ProtectEvidence evidence2 =
                ProtectTestFactory.builder().gene(null).event("Event 2").build();
        evidences.add(evidence2);
        evidences.add(evidence2);
        evidences.add(evidence2);

        assertEquals(2, EventKey.buildUniqueEventSet(evidences).size());
    }
}