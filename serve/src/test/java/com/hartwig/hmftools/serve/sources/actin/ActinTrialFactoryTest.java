package com.hartwig.hmftools.serve.sources.actin;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;

import org.junit.Test;

public class ActinTrialFactoryTest {

    @Test
    public void canExtractActinTrials() {
        ActinEntry entry1 = ActinTestFactory.builder().trial("trial 1").cohort("A").isUsedAsInclusion(true).build();
        ActinTrial trial1 = ActinTrialFactory.toActinTrial(entry1, "some event");
        assertEquals("some event", trial1.sourceEvent());
        assertEquals("trial 1|A", trial1.treatment());
        assertEquals(EvidenceDirection.RESPONSIVE, trial1.direction());

        ActinEntry entry2 = ActinTestFactory.builder().trial("trial 2").isUsedAsInclusion(false).build();
        ActinTrial trial2 = ActinTrialFactory.toActinTrial(entry2, "some other event");
        assertEquals("some other event", trial2.sourceEvent());
        assertEquals("trial 2", trial2.treatment());
        assertEquals(EvidenceDirection.NO_BENEFIT, trial2.direction());
    }
}