package com.hartwig.hmftools.protect.evidence;

import static com.hartwig.hmftools.common.serve.actionability.EvidenceDirection.RESPONSIVE;
import static com.hartwig.hmftools.common.serve.actionability.EvidenceLevel.A;
import static com.hartwig.hmftools.common.serve.actionability.EvidenceLevel.B;
import static com.hartwig.hmftools.common.serve.actionability.EvidenceLevel.C;
import static com.hartwig.hmftools.protect.evidence.ProtectEvidenceItems.highestReportableLevel;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Optional;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidenceItem;
import com.hartwig.hmftools.common.protect.ProtectEvidenceItem;

import org.junit.Test;

public class ProtectEvidenceItemsTest {

    private final ProtectEvidenceItem onLabelResponsiveA = ProtectEvidenceItemTestFactory.createDefault(true, RESPONSIVE, A).build();
    private final ProtectEvidenceItem onLabelResponsiveB = ProtectEvidenceItemTestFactory.createDefault(true, RESPONSIVE, B).build();
    private final ProtectEvidenceItem onLabelResponsiveC = ProtectEvidenceItemTestFactory.createDefault(true, RESPONSIVE, C).build();
    private final ProtectEvidenceItem offLabelResponsiveA = ProtectEvidenceItemTestFactory.createDefault(false, RESPONSIVE, A).build();
    private final ProtectEvidenceItem offLabelResponsiveB = ProtectEvidenceItemTestFactory.createDefault(false, RESPONSIVE, B).build();
    private final ProtectEvidenceItem offLabelResponsiveC = ProtectEvidenceItemTestFactory.createDefault(false, RESPONSIVE, C).build();

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
        final List<ProtectEvidenceItem> evidence = Lists.newArrayList(onLabelResponsiveC, offLabelResponsiveC);
        final Set<ProtectEvidenceItem> victims = Sets.newHashSet(ProtectEvidenceItems.reportHighest(evidence));
        assertEquals(2, victims.size());
        for (ProtectEvidenceItem victim : victims) {
            assertFalse(victim.reported());
        }
    }

    @Test
    public void doNoReportOffLabelAtSameLevelAsOnLabel() {
        final List<ProtectEvidenceItem> evidence = Lists.newArrayList(onLabelResponsiveA, offLabelResponsiveA);
        final Set<ProtectEvidenceItem> victims = Sets.newHashSet(ProtectEvidenceItems.reportHighest(evidence));
        assertEquals(2, victims.size());
        assertTrue(victims.contains(onLabelResponsiveA));
        assertFalse(victims.contains(offLabelResponsiveA));
    }

    @Test
    public void reportHighestOffLabelIfHigherThanOnLabel() {
        final List<ProtectEvidenceItem> evidence = Lists.newArrayList(onLabelResponsiveC, offLabelResponsiveA, offLabelResponsiveB);
        final Set<ProtectEvidenceItem> victims = Sets.newHashSet(ProtectEvidenceItems.reportHighest(evidence));
        assertEquals(3, victims.size());
        assertTrue(victims.contains(offLabelResponsiveA));
        assertFalse(victims.contains(offLabelResponsiveB));
        assertFalse(victims.contains(onLabelResponsiveC));
    }

    @Test
    public void doNotSetReportToTrue() {
        final ProtectEvidenceItem reported = onLabelResponsiveA;
        final ProtectEvidenceItem reportedVictim = ProtectEvidenceItems.reportHighest(Lists.newArrayList(reported)).get(0);
        assertTrue(reportedVictim.reported());

        final ProtectEvidenceItem notReported = ImmutableProtectEvidenceItem.builder().from(reported).reported(false).build();
        final ProtectEvidenceItem notReportedVictim = ProtectEvidenceItems.reportHighest(Lists.newArrayList(notReported)).get(0);
        assertFalse(notReportedVictim.reported());
    }
}
