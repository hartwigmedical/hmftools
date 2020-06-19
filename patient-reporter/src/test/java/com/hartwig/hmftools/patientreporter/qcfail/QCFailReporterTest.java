package com.hartwig.hmftools.patientreporter.qcfail;

import static org.junit.Assert.*;

import com.hartwig.hmftools.patientreporter.cfreport.ForNumber;

import org.junit.Test;

public class QCFailReporterTest {

    @Test
    public void canDetermineForNumber() {
        assertEquals(ForNumber.FOR_082.display(), QCFailReporter.determineForNumber(QCFailReason.TECHNICAL_FAILURE));
        assertEquals(ForNumber.FOR_083.display(), QCFailReporter.determineForNumber(QCFailReason.SUFFICIENT_TCP_QC_FAILURE));
        assertEquals(ForNumber.FOR_100.display(), QCFailReporter.determineForNumber(QCFailReason.INSUFFICIENT_TCP_SHALLOW_WGS));
        assertEquals(ForNumber.FOR_100.display(), QCFailReporter.determineForNumber(QCFailReason.INSUFFICIENT_TCP_DEEP_WGS));
        assertEquals(ForNumber.FOR_102.display(), QCFailReporter.determineForNumber(QCFailReason.INSUFFICIENT_DNA));
        assertEquals(ForNumber.FOR_UNDEFINED.display(), QCFailReporter.determineForNumber(QCFailReason.UNDEFINED));
        assertNotEquals(ForNumber.FOR_100.display(), QCFailReporter.determineForNumber(QCFailReason.TECHNICAL_FAILURE));
        assertNotEquals(ForNumber.FOR_080.display(), QCFailReporter.determineForNumber(QCFailReason.INSUFFICIENT_DNA));
    }
}