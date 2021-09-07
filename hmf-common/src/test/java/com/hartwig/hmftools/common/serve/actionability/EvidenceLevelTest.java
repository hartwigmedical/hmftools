package com.hartwig.hmftools.common.serve.actionability;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class EvidenceLevelTest {

    @Test
    public void canConvertFromDisplayString() {
        assertEquals(EvidenceLevel.A, EvidenceLevel.fromString("A"));

        assertNull(EvidenceLevel.fromString("don't know what this is"));
    }

    @Test
    public void canCompareEvidenceLevels() {
        assertTrue(EvidenceLevel.A.isHigher(EvidenceLevel.B));
        assertFalse(EvidenceLevel.D.isHigher(EvidenceLevel.C));
        assertFalse(EvidenceLevel.B.isHigher(EvidenceLevel.B));
    }
}