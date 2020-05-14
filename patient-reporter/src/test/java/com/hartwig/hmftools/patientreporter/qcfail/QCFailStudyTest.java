package com.hartwig.hmftools.patientreporter.qcfail;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class QCFailStudyTest {

    @Test
    public void canDeriveStudyFromSampleId() {
        assertEquals(QCFailStudy.CPCT, QCFailStudy.fromSampleId("CPCT02010500T"));
        assertEquals(QCFailStudy.DRUP, QCFailStudy.fromSampleId("DRUP01010400TII"));
        assertEquals(QCFailStudy.WIDE, QCFailStudy.fromSampleId("WIDE01T"));
    }
}