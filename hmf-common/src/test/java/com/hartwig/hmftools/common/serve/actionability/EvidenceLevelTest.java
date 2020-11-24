package com.hartwig.hmftools.common.serve.actionability;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import org.junit.Test;

public class EvidenceLevelTest {

    @Test
    public void canConvertFromDisplayString() {
        assertEquals(EvidenceLevel.A, EvidenceLevel.fromString("A"));

        assertNull(EvidenceLevel.fromString("don't know what this is"));
    }
}