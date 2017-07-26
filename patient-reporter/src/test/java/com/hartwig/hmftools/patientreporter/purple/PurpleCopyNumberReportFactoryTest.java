package com.hartwig.hmftools.patientreporter.purple;

import static com.hartwig.hmftools.patientreporter.purple.PurpleCopyNumberReportFactory.ABS_GAIN;
import static com.hartwig.hmftools.patientreporter.purple.PurpleCopyNumberReportFactory.ABS_LOSS;
import static com.hartwig.hmftools.patientreporter.purple.PurpleCopyNumberReportFactory.REL_GAIN;
import static com.hartwig.hmftools.patientreporter.purple.PurpleCopyNumberReportFactory.type;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReportType;

import org.junit.Test;

public class PurpleCopyNumberReportFactoryTest {

    @Test
    public void testLoss() {
        assertEquals(CopyNumberReportType.LOSS, type(0, ABS_LOSS));
        assertEquals(CopyNumberReportType.LOSS, type(1, ABS_LOSS));
        assertEquals(CopyNumberReportType.LOSS, type(2, ABS_LOSS));
        assertEquals(CopyNumberReportType.LOSS, type(2, ABS_LOSS - 0.1));
        assertEquals(CopyNumberReportType.NEUTRAL, type(2, ABS_LOSS + 0.1));
    }

    @Test
    public void testAbsGain() {
        assertEquals(CopyNumberReportType.GAIN, type(0, ABS_GAIN));
        assertEquals(CopyNumberReportType.GAIN, type(1, ABS_GAIN));
        assertEquals(CopyNumberReportType.GAIN, type(2, ABS_GAIN));
        assertEquals(CopyNumberReportType.GAIN, type(ABS_GAIN, ABS_GAIN));

        assertEquals(CopyNumberReportType.GAIN, type(ABS_GAIN, ABS_GAIN + 0.1));
        assertEquals(CopyNumberReportType.NEUTRAL, type(ABS_GAIN, ABS_GAIN - 0.1));
    }

    @Test
    public void testRelGain() {
        assertEquals(CopyNumberReportType.GAIN, type(0, REL_GAIN));
        assertEquals(CopyNumberReportType.GAIN, type(1, REL_GAIN));
        assertEquals(CopyNumberReportType.GAIN, type(1, REL_GAIN + 0.1));
        assertEquals(CopyNumberReportType.NEUTRAL, type(1, REL_GAIN - 0.1));
        assertEquals(CopyNumberReportType.NEUTRAL, type(1.01, REL_GAIN));
    }
}
