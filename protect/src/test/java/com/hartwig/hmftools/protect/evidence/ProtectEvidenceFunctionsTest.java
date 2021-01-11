package com.hartwig.hmftools.protect.evidence;

import static com.hartwig.hmftools.common.serve.actionability.EvidenceDirection.RESPONSIVE;
import static com.hartwig.hmftools.common.serve.actionability.EvidenceLevel.A;
import static com.hartwig.hmftools.common.serve.actionability.EvidenceLevel.B;
import static com.hartwig.hmftools.common.serve.actionability.EvidenceLevel.C;
import static com.hartwig.hmftools.protect.evidence.ProtectEvidenceFunctions.highestReportableLevel;

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

import org.junit.Test;

public class ProtectEvidenceFunctionsTest {

    private final ProtectEvidence onLabelResponsiveA = ProtectEvidenceTestFactory.createDefault(true, RESPONSIVE, A).build();
    private final ProtectEvidence onLabelResponsiveB = ProtectEvidenceTestFactory.createDefault(true, RESPONSIVE, B).build();
    private final ProtectEvidence onLabelResponsiveC = ProtectEvidenceTestFactory.createDefault(true, RESPONSIVE, C).build();
    private final ProtectEvidence offLabelResponsiveA = ProtectEvidenceTestFactory.createDefault(false, RESPONSIVE, A).build();
    private final ProtectEvidence offLabelResponsiveB = ProtectEvidenceTestFactory.createDefault(false, RESPONSIVE, B).build();
    private final ProtectEvidence offLabelResponsiveC = ProtectEvidenceTestFactory.createDefault(false, RESPONSIVE, C).build();

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
    public void neverReportC() {
        final List<ProtectEvidence> evidence = Lists.newArrayList(onLabelResponsiveC, offLabelResponsiveC);
        final Set<ProtectEvidence> victims = Sets.newHashSet(ProtectEvidenceFunctions.reportHighest(evidence));
        assertEquals(2, victims.size());
        for (ProtectEvidence victim : victims) {
            assertFalse(victim.reported());
        }
    }

    @Test
    public void doNoReportOffLabelAtSameLevelAsOnLabel() {
        final List<ProtectEvidence> evidence = Lists.newArrayList(onLabelResponsiveA, offLabelResponsiveA);
        final Set<ProtectEvidence> victims = Sets.newHashSet(ProtectEvidenceFunctions.reportHighest(evidence));
        assertEquals(2, victims.size());
        assertTrue(victims.contains(onLabelResponsiveA));
        assertFalse(victims.contains(offLabelResponsiveA));
    }

    @Test
    public void reportHighestOffLabelIfHigherThanOnLabel() {
        final List<ProtectEvidence> evidence = Lists.newArrayList(onLabelResponsiveC, offLabelResponsiveA, offLabelResponsiveB);
        final Set<ProtectEvidence> victims = Sets.newHashSet(ProtectEvidenceFunctions.reportHighest(evidence));
        assertEquals(3, victims.size());
        assertTrue(victims.contains(offLabelResponsiveA));
        assertFalse(victims.contains(offLabelResponsiveB));
        assertFalse(victims.contains(onLabelResponsiveC));
    }

    @Test
    public void doNotSetReportToTrue() {
        final ProtectEvidence reported = onLabelResponsiveA;
        final ProtectEvidence reportedVictim = ProtectEvidenceFunctions.reportHighest(Lists.newArrayList(reported)).get(0);
        assertTrue(reportedVictim.reported());

        final ProtectEvidence notReported = ImmutableProtectEvidence.builder().from(reported).reported(false).build();
        final ProtectEvidence notReportedVictim = ProtectEvidenceFunctions.reportHighest(Lists.newArrayList(notReported)).get(0);
        assertFalse(notReportedVictim.reported());
    }
}
