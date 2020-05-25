package com.hartwig.hmftools.patientreporter.actionability;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestFactory.createTestEvidenceBuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ActionabilitySource;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.actionability.EvidenceScope;
import com.hartwig.hmftools.common.actionability.ImmutableEvidenceItem;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Ignore;
import org.junit.Test;

public class ReportableEvidenceItemFactoryTest {

    @Test
    public void higherLevelWorks() {
        EvidenceItem item1 = createTestEvidenceBuilder().level(EvidenceLevel.LEVEL_A).build();
        EvidenceItem item2 = createTestEvidenceBuilder().level(EvidenceLevel.LEVEL_B).build();

        assertTrue(ReportableEvidenceItemFactory.hasHigherEvidence(item1, item2));
        assertFalse(ReportableEvidenceItemFactory.hasHigherEvidence(item2, item1));
        assertFalse(ReportableEvidenceItemFactory.hasHigherEvidence(item1, item1));
    }

    @Test
    public void canSelectHighest() {
        EvidenceItem item1 = createTestEvidenceBuilder().level(EvidenceLevel.LEVEL_A).isOnLabel(true).build();
        EvidenceItem item2 = createTestEvidenceBuilder().level(EvidenceLevel.LEVEL_B).isOnLabel(false).build();
        EvidenceItem item3 = createTestEvidenceBuilder().level(EvidenceLevel.LEVEL_C).isOnLabel(false).build();

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
    public void hasKnownLevel() {
        EvidenceItem item1 = createTestEvidenceBuilder().level(EvidenceLevel.LEVEL_A).isOnLabel(true).build();
        EvidenceItem item2 = createTestEvidenceBuilder().level(EvidenceLevel.LEVEL_B).isOnLabel(true).build();
        EvidenceItem item3 = createTestEvidenceBuilder().level(EvidenceLevel.LEVEL_C).isOnLabel(false).build();
        EvidenceItem item4 = createTestEvidenceBuilder().level(EvidenceLevel.LEVEL_D).isOnLabel(false).build();
        EvidenceItem item5 = createTestEvidenceBuilder().level(EvidenceLevel.LEVEL_E).isOnLabel(false).build();

        assertTrue(ReportableEvidenceItemFactory.hasReportableEvidenceLevel(item1));
        assertTrue(ReportableEvidenceItemFactory.hasReportableEvidenceLevel(item2));
        assertFalse(ReportableEvidenceItemFactory.hasReportableEvidenceLevel(item3));
        assertFalse(ReportableEvidenceItemFactory.hasReportableEvidenceLevel(item4));
        assertFalse(ReportableEvidenceItemFactory.hasReportableEvidenceLevel(item5));
    }

    @Test
    public void canFilterEvidenceItemsWithFiltering() {
        List<EvidenceItem> evidenceItems = Lists.newArrayList();

        evidenceItems.add(testEvidenceBuilder().event("TP53 p.Val600Glu").drug("Tamoxifen").build());
        evidenceItems.add(testEvidenceBuilder().event("BRAF p.Val600Glu").drug("Cobimetinib + Vemurafenib").build());

        assertEquals(1, ReportableEvidenceItemFactory.extractAllReportableEvidenceItems(evidenceItems).size());
    }

    @Test
    @Ignore
    public void canFilterEvidenceItemsWithoutFiltering() {
        List<EvidenceItem> evidenceItems = Lists.newArrayList();

        evidenceItems.add(testEvidenceBuilder().event("TP53 Deletion").drug("Nivolumab").build());
        evidenceItems.add(testEvidenceBuilder().event("BRAF p.Val600Glu").drug("Tamoxifen").build());

        assertEquals(2, ReportableEvidenceItemFactory.extractAllReportableEvidenceItems(evidenceItems).size());
    }

    @NotNull
    private static ImmutableEvidenceItem.Builder testEvidenceBuilder() {
        return ImmutableEvidenceItem.builder()
                .isOnLabel(true)
                .drugsType(Strings.EMPTY)
                .cancerType(Strings.EMPTY)
                .response(Strings.EMPTY)
                .level(EvidenceLevel.LEVEL_A)
                .reference(Strings.EMPTY)
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.SPECIFIC);
    }
}