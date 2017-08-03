package com.hartwig.hmftools.patientreporter.algo;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import org.junit.Test;

public class NotSequenceableStudyTest {

    @Test
    public void canDeriveFromSampleName() {
        assertEquals(NotSequenceableStudy.CPCT, NotSequenceableStudy.fromSample("CPCT02010500T"));
        assertEquals(NotSequenceableStudy.DRUP, NotSequenceableStudy.fromSample("DRUP01010400TII"));
        assertNull(NotSequenceableStudy.fromSample("PMC010001"));
    }
}