package com.hartwig.hmftools.patientreporter.actionability;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.ActionabilitySource;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.actionability.EvidenceScope;
import com.hartwig.hmftools.common.actionability.ImmutableEvidenceItem;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.patientreporter.PatientReporterTestFactory;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportableEvidenceItemFactoryTest {

    @Test
    public void reportableFactoryWorksForTrivialCase() {
        ProtectEvidence item1 = testProtectEvidenceBuilder().genomicEvent("A")
                .treatment("A")
                .level(com.hartwig.hmftools.common.serve.actionability.EvidenceLevel.A)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CIVIC))
                .build();
        ProtectEvidence item2 = testProtectEvidenceBuilder().genomicEvent("A")
                .treatment("A")
                .level(com.hartwig.hmftools.common.serve.actionability.EvidenceLevel.A)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CGI))
                .build();
        ProtectEvidence item3 = testProtectEvidenceBuilder().genomicEvent("B")
                .treatment("B")
                .level(com.hartwig.hmftools.common.serve.actionability.EvidenceLevel.A)
                .sources(Sets.newHashSet(Knowledgebase.ICLUSION))
                .build();
        ProtectEvidence item4 = testProtectEvidenceBuilder().genomicEvent("C")
                .treatment("C")
                .level(com.hartwig.hmftools.common.serve.actionability.EvidenceLevel.C)
                .sources(Sets.newHashSet(Knowledgebase.VICC_CGI))
                .build();

        List<ProtectEvidence> nonTrials = ReportableEvidenceItemFactory.extractNonTrials(Lists.newArrayList(item1, item2, item3, item4));
        assertTrue(nonTrials.contains(item1));
        assertTrue(nonTrials.contains(item2));
      //  assertFalse(nonTrials.contains(item3));
        assertTrue(nonTrials.contains(item4));
    }

    @Test
    public void higherLevelWorks() {
        EvidenceItem item1 = PatientReporterTestFactory.createTestEvidenceBuilder().level(EvidenceLevel.LEVEL_A).build();
        EvidenceItem item2 = PatientReporterTestFactory.createTestEvidenceBuilder().level(EvidenceLevel.LEVEL_B).build();

        assertTrue(ReportableEvidenceItemFactory.hasHigherEvidence(item1, item2));
        assertFalse(ReportableEvidenceItemFactory.hasHigherEvidence(item2, item1));
        assertFalse(ReportableEvidenceItemFactory.hasHigherEvidence(item1, item1));
    }

    @Test
    public void canSelectHighest() {
        EvidenceItem item1 = PatientReporterTestFactory.createTestEvidenceBuilder().level(EvidenceLevel.LEVEL_A).isOnLabel(true).build();
        EvidenceItem item2 = PatientReporterTestFactory.createTestEvidenceBuilder().level(EvidenceLevel.LEVEL_B).isOnLabel(false).build();
        EvidenceItem item3 = PatientReporterTestFactory.createTestEvidenceBuilder().level(EvidenceLevel.LEVEL_C).isOnLabel(false).build();

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
    public void canFilterCorrectlyOnEvidenceLevel() {
        EvidenceItem item1 = PatientReporterTestFactory.createTestEvidenceBuilder().level(EvidenceLevel.LEVEL_A).build();
        EvidenceItem item2 = PatientReporterTestFactory.createTestEvidenceBuilder().level(EvidenceLevel.LEVEL_B).build();
        EvidenceItem item3 = PatientReporterTestFactory.createTestEvidenceBuilder().level(EvidenceLevel.LEVEL_C).build();
        EvidenceItem item4 = PatientReporterTestFactory.createTestEvidenceBuilder().level(EvidenceLevel.LEVEL_D).build();
        EvidenceItem item5 = PatientReporterTestFactory.createTestEvidenceBuilder().level(EvidenceLevel.LEVEL_E).build();

        assertTrue(ReportableEvidenceItemFactory.hasReportableEvidenceLevel(item1));
        assertTrue(ReportableEvidenceItemFactory.hasReportableEvidenceLevel(item2));
        assertFalse(ReportableEvidenceItemFactory.hasReportableEvidenceLevel(item3));
        assertFalse(ReportableEvidenceItemFactory.hasReportableEvidenceLevel(item4));
        assertFalse(ReportableEvidenceItemFactory.hasReportableEvidenceLevel(item5));
    }

    @Test
    public void canFilterBlacklistedEvidence() {
        List<EvidenceItem> evidenceItems = Lists.newArrayList();

        evidenceItems.add(testEvidenceBuilder().event("TP53 Deletion").drug("Tamoxifen").build());
        evidenceItems.add(testEvidenceBuilder().event("BRAF p.Val600Glu").drug("Cobimetinib + Vemurafenib").build());
        evidenceItems.add(testEvidenceBuilder().event("TP53 Deletion").drug("Nivolumab").build());
        evidenceItems.add(testEvidenceBuilder().event("BRAF p.Val600Glu").drug("Tamoxifen").build());

        assertEquals(3, ReportableEvidenceItemFactory.filterBlacklistedEvidence(evidenceItems).size());
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

    @NotNull
    private static ImmutableProtectEvidence.Builder testProtectEvidenceBuilder() {
        return ImmutableProtectEvidence.builder()
                .germline(false)
                .reported(true)
                .onLabel(true)
                .direction(EvidenceDirection.RESPONSIVE)
                .urls(Sets.newHashSet("iclusion"));
    }
}