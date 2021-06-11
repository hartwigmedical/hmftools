package com.hartwig.hmftools.protect;

import static com.hartwig.hmftools.common.serve.actionability.EvidenceDirection.PREDICTED_RESPONSIVE;
import static com.hartwig.hmftools.common.serve.actionability.EvidenceDirection.RESPONSIVE;
import static com.hartwig.hmftools.common.serve.actionability.EvidenceLevel.A;
import static com.hartwig.hmftools.common.serve.actionability.EvidenceLevel.B;
import static com.hartwig.hmftools.common.serve.actionability.EvidenceLevel.C;
import static com.hartwig.hmftools.protect.EvidenceReportingFunctions.highestReportableLevel;
import static com.hartwig.hmftools.protect.ProtectTestFactory.testEvidenceBuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class EvidenceReportingFunctionsTest {

    private final ProtectEvidence onLabelResponsiveA = createTestEvidence(true, RESPONSIVE, A);
    private final ProtectEvidence onLabelResponsiveB = createTestEvidence(true, RESPONSIVE, B);
    private final ProtectEvidence onLabelResponsiveC = createTestEvidence(true, RESPONSIVE, C);
    private final ProtectEvidence offLabelResponsiveA = createTestEvidence(false, RESPONSIVE, A);
    private final ProtectEvidence offLabelResponsiveB = createTestEvidence(false, RESPONSIVE, B);
    private final ProtectEvidence offLabelResponsiveC = createTestEvidence(false, RESPONSIVE, C);

    @Test
    public void canDetermineHighestReportableLevel() {
        assertEquals(A, highestReportableLevel(true, Lists.newArrayList(onLabelResponsiveA, onLabelResponsiveB, offLabelResponsiveA)));
        assertEquals(B, highestReportableLevel(true, Lists.newArrayList(onLabelResponsiveB, offLabelResponsiveA)));
        assertEquals(A, highestReportableLevel(false, Lists.newArrayList(onLabelResponsiveA, onLabelResponsiveB, offLabelResponsiveA)));
        assertEquals(A, highestReportableLevel(false, Lists.newArrayList(offLabelResponsiveA)));
    }

    @Test
    public void respectMaxReportingLevel() {
        List<ProtectEvidence> hartwigEvidences = Lists.newArrayList(ImmutableProtectEvidence.builder()
                        .from(onLabelResponsiveB)
                        .addSources(Knowledgebase.HARTWIG_CURATED)
                        .build(),
                ImmutableProtectEvidence.builder()
                        .from(onLabelResponsiveC)
                        .direction(RESPONSIVE)
                        .addSources(Knowledgebase.HARTWIG_CURATED)
                        .build(),
                ImmutableProtectEvidence.builder()
                        .from(offLabelResponsiveC)
                        .direction(PREDICTED_RESPONSIVE)
                        .addSources(Knowledgebase.HARTWIG_CURATED)
                        .build());
        List<ProtectEvidence> hartwigFiltered = EvidenceReportingFunctions.applyReportingAlgo(hartwigEvidences);
        assertEquals(3, hartwigFiltered.size());
        assertEquals(0, hartwigFiltered.stream().filter(x -> x.reported()).count());

        List<ProtectEvidence> viccEvidences = Lists.newArrayList(ImmutableProtectEvidence.builder()
                        .from(onLabelResponsiveB)
                        .addSources(Knowledgebase.VICC_CGI)
                        .build(),
                ImmutableProtectEvidence.builder()
                        .from(onLabelResponsiveC)
                        .direction(RESPONSIVE)
                        .addSources(Knowledgebase.VICC_CGI)
                        .build(),
                ImmutableProtectEvidence.builder()
                        .from(offLabelResponsiveC)
                        .direction(PREDICTED_RESPONSIVE)
                        .addSources(Knowledgebase.VICC_CGI)
                        .build());
        List<ProtectEvidence> viccFiltered = EvidenceReportingFunctions.applyReportingAlgo(viccEvidences);
        assertEquals(3, viccFiltered.size());
        assertEquals(1, viccFiltered.stream().filter(x -> x.reported()).count());

        List<ProtectEvidence> ckbEvidences = Lists.newArrayList(ImmutableProtectEvidence.builder()
                        .from(onLabelResponsiveB)
                        .addSources(Knowledgebase.CKB)
                        .build(),
                ImmutableProtectEvidence.builder()
                        .from(onLabelResponsiveC)
                        .direction(RESPONSIVE)
                        .addSources(Knowledgebase.CKB)
                        .build(),
                ImmutableProtectEvidence.builder()
                        .from(offLabelResponsiveC)
                        .direction(PREDICTED_RESPONSIVE)
                        .addSources(Knowledgebase.CKB)
                        .build());
        List<ProtectEvidence> ckbFiltered = EvidenceReportingFunctions.applyReportingAlgo(ckbEvidences);
        assertEquals(3, ckbFiltered.size());
        assertEquals(2, ckbFiltered.stream().filter(x -> x.reported()).count());

    }

    @Test
    public void doNoReportOffLabelAtSameLevelAsOnLabel() {
        List<ProtectEvidence> evidence = Lists.newArrayList(onLabelResponsiveA, offLabelResponsiveA);
        List<ProtectEvidence> filtered = EvidenceReportingFunctions.applyReportingAlgo(evidence);
        assertEquals(2, filtered.size());
        assertTrue(filtered.contains(onLabelResponsiveA));
        assertFalse(filtered.contains(offLabelResponsiveA));
    }

    @Test
    public void reportHighestOffLabelIfHigherThanOnLabel() {
        List<ProtectEvidence> evidence = Lists.newArrayList(onLabelResponsiveC, offLabelResponsiveA, offLabelResponsiveB);
        List<ProtectEvidence> filtered = EvidenceReportingFunctions.applyReportingAlgo(evidence);
        assertEquals(3, filtered.size());
        assertTrue(filtered.contains(offLabelResponsiveA));
        assertFalse(filtered.contains(offLabelResponsiveB));
        assertFalse(filtered.contains(onLabelResponsiveC));
    }

    @Test
    public void neverSetReportToTrue() {
        ProtectEvidence reported = onLabelResponsiveA;
        ProtectEvidence reportedFiltered = EvidenceReportingFunctions.applyReportingAlgo(Lists.newArrayList(reported)).get(0);
        assertTrue(reportedFiltered.reported());

        ProtectEvidence notReported = ImmutableProtectEvidence.builder().from(reported).reported(false).build();
        ProtectEvidence notReportedFiltered = EvidenceReportingFunctions.applyReportingAlgo(Lists.newArrayList(notReported)).get(0);
        assertFalse(notReportedFiltered.reported());

        ProtectEvidence evidence1 =
                ImmutableProtectEvidence.builder().from(createTestEvidence(true, RESPONSIVE, A)).reported(false).build();
        ProtectEvidence evidence2 = createTestEvidence(false, RESPONSIVE, A);

        List<ProtectEvidence> filtered = EvidenceReportingFunctions.applyReportingAlgo(Lists.newArrayList(evidence1, evidence2));
        assertTrue(filtered.contains(evidence1));
        assertTrue(filtered.contains(evidence2));
    }

    @Test
    public void doNotReportOffLabelTrials() {
        String event1 = "event1";
        String event2 = "event2";
        String event3 = "event3";
        String event4 = "event4";

        ProtectEvidence evidence1 =
                testEvidenceBuilder().genomicEvent(event1).addSources(Knowledgebase.ICLUSION).reported(true).onLabel(true).build();
        ProtectEvidence evidence2 =
                testEvidenceBuilder().genomicEvent(event2).addSources(Knowledgebase.ICLUSION).reported(true).onLabel(false).build();
        ProtectEvidence evidence3 =
                testEvidenceBuilder().genomicEvent(event3).addSources(Knowledgebase.VICC_CGI).reported(true).onLabel(false).build();
        ProtectEvidence evidence4 = testEvidenceBuilder().genomicEvent(event4)
                .addSources(Knowledgebase.VICC_CGI, Knowledgebase.ICLUSION)
                .reported(true)
                .onLabel(false)
                .build();

        List<ProtectEvidence> evidence =
                EvidenceReportingFunctions.reportOnLabelTrialsOnly(Lists.newArrayList(evidence1, evidence2, evidence3, evidence4));

        assertEquals(4, evidence.size());
        assertTrue(evidence.contains(evidence1));
        assertTrue(evidence.contains(evidence3));
        assertTrue(evidence.contains(evidence4));

        ProtectEvidence convertedEvidence2 = findByEvent(evidence, event2);
        assertFalse(convertedEvidence2.reported());
    }

    @NotNull
    private static ProtectEvidence findByEvent(@NotNull Iterable<ProtectEvidence> evidences, @NotNull String event) {
        for (ProtectEvidence evidence : evidences) {
            if (evidence.genomicEvent().equals(event)) {
                return evidence;
            }
        }

        throw new IllegalStateException("Could not find evidence with genomic event: " + event);
    }

    @NotNull
    private static ProtectEvidence createTestEvidence(boolean onLabel, @NotNull EvidenceDirection direction, @NotNull EvidenceLevel level) {
        return testEvidenceBuilder().onLabel(onLabel).level(level).direction(direction).build();
    }
}
