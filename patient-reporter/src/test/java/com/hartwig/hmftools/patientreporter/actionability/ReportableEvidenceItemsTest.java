package com.hartwig.hmftools.patientreporter.actionability;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ActionabilitySource;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.ImmutableEvidenceItem;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportableEvidenceItemsTest {

    @Test
    public void higherLevelWorks() {
        EvidenceItem item1 = builder().level("A").build();
        EvidenceItem item2 = builder().level("B").build();

        assertTrue(ReportableEvidenceItems.hasHigherEvidence(item1, item2));
        assertFalse(ReportableEvidenceItems.hasHigherEvidence(item2, item1));
        assertFalse(ReportableEvidenceItems.hasHigherEvidence(item1, item1));
    }

    @Test
    public void canSelectHighest() {
        EvidenceItem item1 = builder().level("A").isOnLabel(true).build();
        EvidenceItem item2 = builder().level("B").isOnLabel(false).build();
        EvidenceItem item3 = builder().level("C").isOnLabel(false).build();

        EvidenceItem highestOffLabel = ReportableEvidenceItems.highestOffLabel(Lists.newArrayList(item1, item2, item3));
        assertNotNull(highestOffLabel);
        assertEquals("B", highestOffLabel.level());


        EvidenceItem highestOnLabel = ReportableEvidenceItems.highestOnLabel(Lists.newArrayList(item1, item2, item3));
        assertNotNull(highestOnLabel);
        assertEquals("A", highestOnLabel.level());

        EvidenceItem onLabelOnly = ReportableEvidenceItems.highestOffLabel(Lists.newArrayList(item1));
        assertNull(onLabelOnly);
    }

    @NotNull
    private static ImmutableEvidenceItem.Builder builder() {
        return ImmutableEvidenceItem.builder()
                .event("")
                .isOnLabel(false)
                .response("")
                .level("")
                .drugsType("")
                .drug("")
                .reference("")
                .source(ActionabilitySource.CIVIC);
    }
}