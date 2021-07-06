package com.hartwig.hmftools.common.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;

import org.junit.Test;

public class PurpleDataLoaderTest {

    @Test
    public void canSeparateReportedFromUnreportedGainsLosses() {
        ReportableGainLoss gain = PurpleTestFactory.testReportableGainLoss("gene1", CopyNumberInterpretation.FULL_GAIN);
        ReportableGainLoss loss = PurpleTestFactory.testReportableGainLoss("gene2", CopyNumberInterpretation.FULL_LOSS);

        List<ReportableGainLoss> reportable = Lists.newArrayList(gain);
        List<ReportableGainLoss> all = Lists.newArrayList(gain, loss);

        List<ReportableGainLoss> unreported = PurpleDataLoader.selectUnreportedGainsLosses(all, reportable);
        assertEquals(1, unreported.size());
        assertTrue(unreported.contains(loss));
        assertFalse(unreported.contains(gain));
    }
}