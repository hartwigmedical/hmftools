package com.hartwig.hmftools.orange.conversion;

import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.common.chord.ChordTestFactory;
import com.hartwig.hmftools.common.doid.DoidTestFactory;
import com.hartwig.hmftools.common.metrics.FlagstatTestFactory;
import com.hartwig.hmftools.common.metrics.BamMetricsTestFactory;
import com.hartwig.hmftools.common.peach.PeachTestFactory;

import org.junit.Test;

public class OrangeConversionTest
{
    @Test
    public void shouldConvertTestVersionsOfAllDatamodels()
    {
        assertNotNull(OrangeConversion.convert(FlagstatTestFactory.createMinimalTestFlagstat()));
        assertNotNull(OrangeConversion.convert(BamMetricsTestFactory.createMinimalTestWGSMetrics()));
        assertNotNull(OrangeConversion.convert(DoidTestFactory.createTestDoidNode()));
        assertNotNull(OrangeConversion.convert(ChordTestFactory.createMinimalTestChordAnalysis()));
        assertNotNull(OrangeConversion.convert(PeachTestFactory.builder().build()));
    }
}
