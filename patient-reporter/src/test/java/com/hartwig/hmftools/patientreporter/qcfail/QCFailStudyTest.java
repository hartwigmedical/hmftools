package com.hartwig.hmftools.patientreporter.qcfail;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class QCFailStudyTest {

    @Test
    public void canDeriveStudyFromSampleName() {
        assertEquals(QCFailStudy.CPCT, QCFailStudy.fromSample("CPCT02010500T"));
        assertEquals(QCFailStudy.DRUP, QCFailStudy.fromSample("DRUP01010400TII"));
        assertEquals(QCFailStudy.WIDE, QCFailStudy.fromSample("WIDE01T"));
    }
}