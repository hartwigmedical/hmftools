package com.hartwig.hmftools.common.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.interpretation.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.interpretation.GainLoss;
import com.hartwig.hmftools.common.purple.interpretation.GainLossTestFactory;

import org.junit.Test;

public class PurpleDataLoaderTest {

    @Test
    public void canSeparateReportedFromUnreportedGainsLosses() {
        GainLoss gain = GainLossTestFactory.createGainLoss("gene1", CopyNumberInterpretation.FULL_GAIN);
        GainLoss loss = GainLossTestFactory.createGainLoss("gene2", CopyNumberInterpretation.FULL_LOSS);

        List<GainLoss> reportable = Lists.newArrayList(gain);
        List<GainLoss> all = Lists.newArrayList(gain, loss);

        List<GainLoss> unreported = PurpleDataLoader.selectUnreportedGainsLosses(all, reportable);
        assertEquals(1, unreported.size());
        assertTrue(unreported.contains(loss));
        assertFalse(unreported.contains(gain));
    }
}