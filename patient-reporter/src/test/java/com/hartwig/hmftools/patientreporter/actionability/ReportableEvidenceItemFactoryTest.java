package com.hartwig.hmftools.patientreporter.actionability;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ActionabilitySource;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.actionability.ImmutableEvidenceItem;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportableEvidenceItemFactoryTest {

    @Test
    public void higherLevelWorks() {
        EvidenceItem item1 = builder().level(EvidenceLevel.LEVELA).build();
        EvidenceItem item2 = builder().level(EvidenceLevel.LEVELB).build();

        assertTrue(ReportableEvidenceItemFactory.hasHigherEvidence(item1, item2));
        assertFalse(ReportableEvidenceItemFactory.hasHigherEvidence(item2, item1));
        assertFalse(ReportableEvidenceItemFactory.hasHigherEvidence(item1, item1));
    }

    @Test
    public void canSelectHighest() {
        EvidenceItem item1 = builder().level(EvidenceLevel.LEVELA).isOnLabel(true).build();
        EvidenceItem item2 = builder().level(EvidenceLevel.LEVELB).isOnLabel(false).build();
        EvidenceItem item3 = builder().level(EvidenceLevel.LEVELC).isOnLabel(false).build();

        EvidenceItem highestOffLabel = ReportableEvidenceItemFactory.highestOffLabel(Lists.newArrayList(item1, item2, item3));
        assertNotNull(highestOffLabel);
        assertEquals(EvidenceLevel.LEVELB, highestOffLabel.level());


        EvidenceItem highestOnLabel = ReportableEvidenceItemFactory.highestOnLabel(Lists.newArrayList(item1, item2, item3));
        assertNotNull(highestOnLabel);
        assertEquals(EvidenceLevel.LEVELA, highestOnLabel.level());

        EvidenceItem onLabelOnly = ReportableEvidenceItemFactory.highestOffLabel(Lists.newArrayList(item1));
        assertNull(onLabelOnly);
    }

    @NotNull
    private static ImmutableEvidenceItem.Builder builder() {
        return ImmutableEvidenceItem.builder()
                .event(Strings.EMPTY)
                .source(ActionabilitySource.CIVIC)
                .reference(Strings.EMPTY)
                .drug(Strings.EMPTY)
                .drugsType(Strings.EMPTY)
                .level(EvidenceLevel.LEVELA)
                .response(Strings.EMPTY)
                .isOnLabel(false);
    }
}