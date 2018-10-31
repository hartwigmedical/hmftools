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
        EvidenceItem item1 = builder().level(EvidenceLevel.LEVEL_A).build();
        EvidenceItem item2 = builder().level(EvidenceLevel.LEVEL_B).build();

        assertTrue(ReportableEvidenceItemFactory.hasHigherEvidence(item1, item2));
        assertFalse(ReportableEvidenceItemFactory.hasHigherEvidence(item2, item1));
        assertFalse(ReportableEvidenceItemFactory.hasHigherEvidence(item1, item1));
    }

    @Test
    public void canSelectHighest() {
        EvidenceItem item1 = builder().level(EvidenceLevel.LEVEL_A).isOnLabel(true).build();
        EvidenceItem item2 = builder().level(EvidenceLevel.LEVEL_B).isOnLabel(false).build();
        EvidenceItem item3 = builder().level(EvidenceLevel.LEVEL_C).isOnLabel(false).build();

        EvidenceItem highestOffLabel = ReportableEvidenceItemFactory.highestOffLabel(Lists.newArrayList(item1, item2, item3));
        assertNotNull(highestOffLabel);
        assertEquals(EvidenceLevel.LEVEL_B, highestOffLabel.level());


        EvidenceItem highestOnLabel = ReportableEvidenceItemFactory.highestOnLabel(Lists.newArrayList(item1, item2, item3));
        assertNotNull(highestOnLabel);
        assertEquals(EvidenceLevel.LEVEL_A, highestOnLabel.level());

        EvidenceItem onLabelOnly = ReportableEvidenceItemFactory.highestOffLabel(Lists.newArrayList(item1));
        assertNull(onLabelOnly);
    }

    @Test
    public void isNotTrail() {
        EvidenceItem item1 = builder().source(ActionabilitySource.ICLUSION).build();
        EvidenceItem item2 = builder().source(ActionabilitySource.ONCOKB).build();
        EvidenceItem item3 = builder().source(ActionabilitySource.CGI).build();
        EvidenceItem item4 = builder().source(ActionabilitySource.CIVIC).build();

        assertFalse(ReportableEvidenceItemFactory.filterNonTrial(item1));
        assertTrue(ReportableEvidenceItemFactory.filterNonTrial(item2));
        assertTrue(ReportableEvidenceItemFactory.filterNonTrial(item3));
        assertTrue(ReportableEvidenceItemFactory.filterNonTrial(item4));
    }

    @Test
    public void hasKnownLevel() {
        EvidenceItem item1 = builder().level(EvidenceLevel.LEVEL_A).isOnLabel(true).build();
        EvidenceItem item2 = builder().level(EvidenceLevel.LEVEL_B).isOnLabel(true).build();
        EvidenceItem item3 = builder().level(EvidenceLevel.LEVEL_C).isOnLabel(false).build();
        EvidenceItem item4 = builder().level(EvidenceLevel.LEVEL_D).isOnLabel(false).build();
        EvidenceItem item5 = builder().level(EvidenceLevel.LEVEL_E).isOnLabel(false).build();

        assertTrue(ReportableEvidenceItemFactory.selectLevelsAandB(item1));
        assertTrue(ReportableEvidenceItemFactory.selectLevelsAandB(item2));
        assertFalse(ReportableEvidenceItemFactory.selectLevelsAandB(item3));
        assertFalse(ReportableEvidenceItemFactory.selectLevelsAandB(item4));
        assertFalse(ReportableEvidenceItemFactory.selectLevelsAandB(item5));
    }

    @NotNull
    private static ImmutableEvidenceItem.Builder builder() {
        return ImmutableEvidenceItem.builder()
                .event(Strings.EMPTY)
                .source(ActionabilitySource.CIVIC)
                .reference(Strings.EMPTY)
                .drug(Strings.EMPTY)
                .drugsType(Strings.EMPTY)
                .level(EvidenceLevel.LEVEL_A)
                .response(Strings.EMPTY)
                .isOnLabel(false);
    }
}