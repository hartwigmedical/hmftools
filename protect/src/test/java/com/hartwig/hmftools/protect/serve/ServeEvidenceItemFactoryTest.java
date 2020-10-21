package com.hartwig.hmftools.protect.serve;

import static com.hartwig.hmftools.protect.serve.ServeEvidenceItemFactory.minLevel;
import static com.hartwig.hmftools.protect.serve.ServeEvidenceItemFactory.report;
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

public class ServeEvidenceItemFactoryTest {

    @Test
    public void testMinLevel() {
        ServeEvidenceItem onLabelResponsiveA = ServeEvidenceItemTest.createDefault(true, RESPONSIVE, A).build();
        ServeEvidenceItem onLabelResponsiveB = ServeEvidenceItemTest.createDefault(true, RESPONSIVE, B).build();
        ServeEvidenceItem offLabelResponsiveA = ServeEvidenceItemTest.createDefault(false, RESPONSIVE, A).build();

        assertEquals(A, minLevel(RESPONSIVE, Lists.newArrayList(onLabelResponsiveA, onLabelResponsiveB, offLabelResponsiveA)));
        assertEquals(B, minLevel(RESPONSIVE, Lists.newArrayList(onLabelResponsiveB, offLabelResponsiveA)));
        assertEquals(A, minLevel(RESPONSIVE, Lists.newArrayList(offLabelResponsiveA)));
        assertEquals(D, minLevel(RESPONSIVE, Lists.newArrayList()));

        assertEquals(D, minLevel(RESISTANT, Lists.newArrayList()));
        assertEquals(D, minLevel(RESISTANT, Lists.newArrayList(onLabelResponsiveA, onLabelResponsiveB, offLabelResponsiveA)));
    }

    @Test
    public void testReport() {
        ServeEvidenceItem responsiveB = ServeEvidenceItemTest.createDefault(true, RESPONSIVE, B).build();
        ServeEvidenceItem resistantC = ServeEvidenceItemTest.createDefault(true, RESISTANT, C).build();

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
