package com.hartwig.hmftools.protect.evidence;

import static com.hartwig.hmftools.protect.evidence.ProtectEvidenceItems.highestLevel;
import static com.hartwig.hmftools.protect.evidence.ProtectEvidenceItems.report;
import static com.hartwig.hmftools.serve.actionability.EvidenceDirection.RESISTANT;
import static com.hartwig.hmftools.serve.actionability.EvidenceDirection.RESPONSIVE;
import static com.hartwig.hmftools.serve.actionability.EvidenceLevel.A;
import static com.hartwig.hmftools.serve.actionability.EvidenceLevel.B;
import static com.hartwig.hmftools.serve.actionability.EvidenceLevel.C;
import static com.hartwig.hmftools.serve.actionability.EvidenceLevel.D;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;

import org.junit.Test;

public class ProtectEvidenceItemsTest {

    @Test
    public void testHighestLevel() {
        ProtectEvidenceItem onLabelResponsiveA = ProtectEvidenceItemTest.createDefault(true, RESPONSIVE, A).build();
        ProtectEvidenceItem onLabelResponsiveB = ProtectEvidenceItemTest.createDefault(true, RESPONSIVE, B).build();
        ProtectEvidenceItem offLabelResponsiveA = ProtectEvidenceItemTest.createDefault(false, RESPONSIVE, A).build();

        assertEquals(A, highestLevel(RESPONSIVE, Lists.newArrayList(onLabelResponsiveA, onLabelResponsiveB, offLabelResponsiveA)));
        assertEquals(B, highestLevel(RESPONSIVE, Lists.newArrayList(onLabelResponsiveB, offLabelResponsiveA)));
        assertEquals(A, highestLevel(RESPONSIVE, Lists.newArrayList(offLabelResponsiveA)));
        assertEquals(D, highestLevel(RESPONSIVE, Lists.newArrayList()));

        assertEquals(D, highestLevel(RESISTANT, Lists.newArrayList()));
        assertEquals(D, highestLevel(RESISTANT, Lists.newArrayList(onLabelResponsiveA, onLabelResponsiveB, offLabelResponsiveA)));
    }

    @Test
    public void testReport() {
        ProtectEvidenceItem responsiveB = ProtectEvidenceItemTest.createDefault(true, RESPONSIVE, B).build();
        ProtectEvidenceItem resistantC = ProtectEvidenceItemTest.createDefault(true, RESISTANT, C).build();

        assertFalse(report(responsiveB, A, A));
        assertTrue(report(responsiveB, B, A));
        assertTrue(report(responsiveB, C, A));
        assertTrue(report(responsiveB, D, A));

        assertFalse(report(resistantC, A, A));
        assertFalse(report(resistantC, A, B));
        assertTrue(report(resistantC, A, C));
        assertTrue(report(resistantC, A, D));
    }

}
