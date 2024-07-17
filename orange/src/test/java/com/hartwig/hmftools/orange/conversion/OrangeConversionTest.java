package com.hartwig.hmftools.orange.conversion;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.common.chord.ChordTestFactory;
import com.hartwig.hmftools.common.doid.DoidTestFactory;
import com.hartwig.hmftools.common.flagstat.FlagstatTestFactory;
import com.hartwig.hmftools.common.lilac.LilacTestFactory;
import com.hartwig.hmftools.common.metrics.WGSMetricsTestFactory;
import com.hartwig.hmftools.common.peach.PeachTestFactory;
import com.hartwig.hmftools.common.virus.VirusTestFactory;
import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.datamodel.virus.VirusInterpretation;

import org.junit.Test;

public class OrangeConversionTest
{
    @Test
    public void shouldConvertTestVersionsOfAllDatamodels()
    {
        assertNotNull(OrangeConversion.convert(FlagstatTestFactory.createMinimalTestFlagstat()));
        assertNotNull(OrangeConversion.convert(WGSMetricsTestFactory.createMinimalTestWGSMetrics()));
        assertNotNull(OrangeConversion.convert(DoidTestFactory.createTestDoidNode()));
        assertNotNull(OrangeConversion.convert(LilacTestFactory.createEmptyData(), true, true));
        assertNotNull(OrangeConversion.convert(LilacTestFactory.alleleBuilder().build(), true, true));
        assertNotNull(OrangeConversion.convert(VirusTestFactory.annotatedVirusBuilder().build()));
        assertNotNull(OrangeConversion.convert(ChordTestFactory.createMinimalTestChordAnalysis()));
        assertNotNull(OrangeConversion.convert(PeachTestFactory.builder().build()));
    }

    @Test
    public void shouldNullLilacFragmentsIfUnavailable()
    {
        assertNotNull(OrangeConversion.convert(LilacTestFactory.alleleBuilder().build(), true, true).refFragments());
        assertNotNull(OrangeConversion.convert(LilacTestFactory.alleleBuilder().build(), true, true).rnaFragments());
        assertNull(OrangeConversion.convert(LilacTestFactory.alleleBuilder().build(), false, true).refFragments());
        assertNull(OrangeConversion.convert(LilacTestFactory.alleleBuilder().build(), true, false).rnaFragments());
    }

    @Test
    public void convertsEachVirusConstantProperly()
    {
        for(final VirusType value : VirusType.values())
        {
            VirusInterpretation.valueOf(value.name());
        }
    }
}