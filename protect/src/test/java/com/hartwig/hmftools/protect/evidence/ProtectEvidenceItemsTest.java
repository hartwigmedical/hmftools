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

    private final ProtectEvidenceItem onLabelResponsiveA = ProtectEvidenceItemTest.createDefault(true, RESPONSIVE, A).build();
    private final ProtectEvidenceItem onLabelResponsiveB = ProtectEvidenceItemTest.createDefault(true, RESPONSIVE, B).build();
    private final ProtectEvidenceItem onLabelResponsiveC = ProtectEvidenceItemTest.createDefault(true, RESPONSIVE, C).build();
    private final ProtectEvidenceItem offLabelResponsiveA = ProtectEvidenceItemTest.createDefault(false, RESPONSIVE, A).build();
    private final ProtectEvidenceItem offLabelResponsiveB = ProtectEvidenceItemTest.createDefault(false, RESPONSIVE, B).build();
    private final ProtectEvidenceItem offLabelResponsiveC = ProtectEvidenceItemTest.createDefault(false, RESPONSIVE, C).build();

    @Test
    public void testHighestReportableLevel() {
        assertEquals(Optional.of(A), highestReportableLevel(true, Lists.newArrayList(onLabelResponsiveA, onLabelResponsiveB, offLabelResponsiveA)));
        assertEquals(Optional.of(B), highestReportableLevel(true, Lists.newArrayList(onLabelResponsiveB, offLabelResponsiveA)));
        assertEquals(Optional.of(A), highestReportableLevel(false, Lists.newArrayList(onLabelResponsiveA, onLabelResponsiveB, offLabelResponsiveA)));
        assertEquals(Optional.of(A), highestReportableLevel(false, Lists.newArrayList(offLabelResponsiveA)));
    }

    @Test
    public void testDoNotReportC() {
        final List<ProtectEvidenceItem> evidence = Lists.newArrayList(onLabelResponsiveC, offLabelResponsiveC);
        final Set<ProtectEvidenceItem> victims = Sets.newHashSet(ProtectEvidenceItems.reportHighest(evidence));
        assertEquals(2, victims.size());
        assertFalse(victims.contains(onLabelResponsiveC));
        assertFalse(victims.contains(offLabelResponsiveC));
    }

    @Test
    public void testDoNoReportOffLabelAtSameLevelAsOnLabel() {
        final List<ProtectEvidenceItem> evidence = Lists.newArrayList(onLabelResponsiveA, offLabelResponsiveA);
        final Set<ProtectEvidenceItem> victims = Sets.newHashSet(ProtectEvidenceItems.reportHighest(evidence));
        assertEquals(2, victims.size());
        assertTrue(victims.contains(onLabelResponsiveA));
        assertFalse(victims.contains(offLabelResponsiveA));
    }

    @Test
    public void testReportHighestOffLabelIfHigherThanOnLabel() {
        final List<ProtectEvidenceItem> evidence = Lists.newArrayList(onLabelResponsiveC, offLabelResponsiveA, offLabelResponsiveB);
        final Set<ProtectEvidenceItem> victims = Sets.newHashSet(ProtectEvidenceItems.reportHighest(evidence));
        assertEquals(3, victims.size());
        assertTrue(victims.contains(offLabelResponsiveA));
        assertFalse(victims.contains(offLabelResponsiveB));
        assertFalse(victims.contains(onLabelResponsiveC));
    }

    @Test
    public void testDoNotSetReportToTrue() {
        final ProtectEvidenceItem reported = onLabelResponsiveA;
        final ProtectEvidenceItem reportedVictim = ProtectEvidenceItems.reportHighest(Lists.newArrayList(reported)).get(0);
        assertTrue(reportedVictim.reported());

        final ProtectEvidenceItem notReported = ImmutableProtectEvidenceItem.builder().from(reported).reported(false).build();
        final ProtectEvidenceItem notReportedVictim = ProtectEvidenceItems.reportHighest(Lists.newArrayList(notReported)).get(0);
        assertFalse(notReportedVictim.reported());
    }

}
