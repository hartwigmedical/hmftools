package com.hartwig.hmftools.patientreporter.copynumber;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class CopyNumberReportTest {

    @Test
    public void canResolveType() {
        final CopyNumberReport gain = new CopyNumberReport.Builder().copyNumber(5).build();
        assertEquals(CopyNumberReport.COPY_NUMBER_GAIN, gain.resolveType());

        final CopyNumberReport loss = new CopyNumberReport.Builder().copyNumber(0).build();
        assertEquals(CopyNumberReport.COPY_NUMBER_LOSS, loss.resolveType());

        final CopyNumberReport neutral = new CopyNumberReport.Builder().copyNumber(2).build();
        assertEquals(CopyNumberReport.COPY_NUMBER_NEUTRAL, neutral.resolveType());
    }
}