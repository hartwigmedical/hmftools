package com.hartwig.hmftools.protect;

import static com.hartwig.hmftools.common.serve.actionability.EvidenceDirection.RESPONSIVE;
import static com.hartwig.hmftools.common.serve.actionability.EvidenceLevel.A;
import static com.hartwig.hmftools.common.serve.actionability.EvidenceLevel.B;
import static com.hartwig.hmftools.common.serve.actionability.EvidenceLevel.C;
import static com.hartwig.hmftools.protect.EvidenceReportingFunctions.highestReportableLevel;
import static com.hartwig.hmftools.protect.ProtectTestFactory.createTestBuilder;
import static com.hartwig.hmftools.protect.ProtectTestFactory.createTestEvidence;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Optional;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class EvidenceReportingFunctionsTest {

    private final ProtectEvidence onLabelResponsiveA = createTestEvidence(true, RESPONSIVE, A).build();
    private final ProtectEvidence onLabelResponsiveB = createTestEvidence(true, RESPONSIVE, B).build();
    private final ProtectEvidence onLabelResponsiveC = createTestEvidence(true, RESPONSIVE, C).build();
    private final ProtectEvidence offLabelResponsiveA = createTestEvidence(false, RESPONSIVE, A).build();
    private final ProtectEvidence offLabelResponsiveB = createTestEvidence(false, RESPONSIVE, B).build();
    private final ProtectEvidence offLabelResponsiveC = createTestEvidence(false, RESPONSIVE, C).build();

    @Test
    public void canDetermineHighestReportableLevel() {
        assertEquals(Optional.of(A),
                highestReportableLevel(true, Lists.newArrayList(onLabelResponsiveA, onLabelResponsiveB, offLabelResponsiveA)));
        assertEquals(Optional.of(B), highestReportableLevel(true, Lists.newArrayList(onLabelResponsiveB, offLabelResponsiveA)));
        assertEquals(Optional.of(A),
                highestReportableLevel(false, Lists.newArrayList(onLabelResponsiveA, onLabelResponsiveB, offLabelResponsiveA)));
        assertEquals(Optional.of(A), highestReportableLevel(false, Lists.newArrayList(offLabelResponsiveA)));
    }

    @Test
    public void doNotReportC() {
        List<ProtectEvidence> evidence = Lists.newArrayList(onLabelResponsiveC, offLabelResponsiveC);
        Set<ProtectEvidence> victims = Sets.newHashSet(EvidenceReportingFunctions.reportHighestLevelEvidence(evidence));
        assertEquals(2, victims.size());
        for (ProtectEvidence victim : victims) {
            assertFalse(victim.reported());
        }
    }

    @Test
    public void doNoReportOffLabelAtSameLevelAsOnLabel() {
        List<ProtectEvidence> evidence = Lists.newArrayList(onLabelResponsiveA, offLabelResponsiveA);
        Set<ProtectEvidence> victims = Sets.newHashSet(EvidenceReportingFunctions.reportHighestLevelEvidence(evidence));
        assertEquals(2, victims.size());
        assertTrue(victims.contains(onLabelResponsiveA));
        assertFalse(victims.contains(offLabelResponsiveA));
    }

    @Test
    public void reportHighestOffLabelIfHigherThanOnLabel() {
        List<ProtectEvidence> evidence = Lists.newArrayList(onLabelResponsiveC, offLabelResponsiveA, offLabelResponsiveB);
        Set<ProtectEvidence> victims = Sets.newHashSet(EvidenceReportingFunctions.reportHighestLevelEvidence(evidence));
        assertEquals(3, victims.size());
        assertTrue(victims.contains(offLabelResponsiveA));
        assertFalse(victims.contains(offLabelResponsiveB));
        assertFalse(victims.contains(onLabelResponsiveC));
    }

    @Test
    public void neverSetReportToTrue() {
        ProtectEvidence reported = onLabelResponsiveA;
        ProtectEvidence reportedVictim = EvidenceReportingFunctions.reportHighestLevelEvidence(Lists.newArrayList(reported)).get(0);
        assertTrue(reportedVictim.reported());

        ProtectEvidence notReported = ImmutableProtectEvidence.builder().from(reported).reported(false).build();
        ProtectEvidence notReportedVictim = EvidenceReportingFunctions.reportHighestLevelEvidence(Lists.newArrayList(notReported)).get(0);
        assertFalse(notReportedVictim.reported());
    }

    @Test
    public void doNotReportOffLabelTrials() {
        String event1 = "event1";
        String event2 = "event2";
        String event3 = "event3";
        String event4 = "event4";

        ProtectEvidence evidence1 =
                createTestBuilder().genomicEvent(event1).addSources(Knowledgebase.ICLUSION).reported(true).onLabel(true).build();
        ProtectEvidence evidence2 =
                createTestBuilder().genomicEvent(event2).addSources(Knowledgebase.ICLUSION).reported(true).onLabel(false).build();
        ProtectEvidence evidence3 =
                createTestBuilder().genomicEvent(event3).addSources(Knowledgebase.VICC_CGI).reported(true).onLabel(false).build();
        ProtectEvidence evidence4 = createTestBuilder().genomicEvent(event4)
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
}
