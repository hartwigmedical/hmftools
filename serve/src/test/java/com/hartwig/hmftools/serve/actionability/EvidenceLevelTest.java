package com.hartwig.hmftools.serve.actionability;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import org.junit.Test;

public class EvidenceLevelTest {

    @Test
    public void canConvertFromString() {
        assertEquals(EvidenceLevel.D, EvidenceLevel.fromString("D"));
        assertNull(EvidenceLevel.fromString("XXX"));
    }
}