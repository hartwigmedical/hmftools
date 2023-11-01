package com.hartwig.hmftools.orange.conversion;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.chord.ChordTestFactory;
import com.hartwig.hmftools.common.doid.DoidTestFactory;
import com.hartwig.hmftools.common.flagstat.FlagstatTestFactory;
import com.hartwig.hmftools.common.lilac.LilacTestFactory;
import com.hartwig.hmftools.common.metrics.WGSMetricsTestFactory;
import com.hartwig.hmftools.common.peach.PeachTestFactory;
import com.hartwig.hmftools.common.sigs.SignatureTestFactory;
import com.hartwig.hmftools.common.virus.VirusTestFactory;
import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.datamodel.virus.VirusInterpretation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class OrangeConversionTest
{
    private static final double EPSILON = 0.01;

    @Test
    public void shouldConvertTestVersionsOfAllDatamodels()
    {
        assertNotNull(OrangeConversion.convert(FlagstatTestFactory.createMinimalTestFlagstat()));
        assertNotNull(OrangeConversion.convert(WGSMetricsTestFactory.createMinimalTestWGSMetrics()));
        assertNotNull(OrangeConversion.convert(DoidTestFactory.createTestDoidNode()));
        assertNotNull(OrangeConversion.convert(LilacTestFactory.createEmptyData(), true, true));
        assertNotNull(OrangeConversion.convert(LilacTestFactory.alleleBuilder().build(), true, true));
        assertNotNull(OrangeConversion.convert(VirusTestFactory.createEmptyData()));
        assertNotNull(OrangeConversion.convert(VirusTestFactory.annotatedVirusBuilder().build()));
        assertNotNull(OrangeConversion.convert(ChordTestFactory.createMinimalTestChordAnalysis()));
        assertNotNull(OrangeConversion.convert(PeachTestFactory.builder().build()));
        assertNotNull(OrangeConversion.convert(SignatureTestFactory.builder().build()));
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
    public void shouldConvertAnnotatedVirus()
    {
        assertTrue(OrangeConversion.convert(VirusTestFactory.createEmptyData()).allViruses().isEmpty());
        assertEqualsValue(VirusTestFactory.createMinimalData(), OrangeConversion.convert(VirusTestFactory.createMinimalData()));
        assertEqualsValue(VirusTestFactory.createProperData(), OrangeConversion.convert(VirusTestFactory.createProperData()));
        assertEqualsValue(VirusTestFactory.createHHVInterpretationData(), OrangeConversion.convert(VirusTestFactory.createHHVInterpretationData()));
    }

    @Test
    public void convertsEachVirusConstantProperly()
    {
        for(final VirusType value : VirusType.values())
        {
            VirusInterpretation.valueOf(value.name());
        }
    }

    private static void assertEqualsValue(@NotNull com.hartwig.hmftools.common.virus.VirusInterpreterData input,
            @NotNull VirusInterpreterData converted)
    {
        assertEqualsValue(input.allViruses(), converted.allViruses());
        assertEqualsValue(input.reportableViruses(), converted.reportableViruses());
    }

    private static void assertEqualsValue(@NotNull List<com.hartwig.hmftools.common.virus.AnnotatedVirus> input,
            @NotNull List<VirusInterpreterEntry> converted)
    {
        assertEquals(converted.size(), input.size());
        for(int i = 0; i < input.size(); i++)
        {
            assertEqualsValue(input.get(i), converted.get(i));
        }
    }

    private static void assertEqualsValue(@NotNull com.hartwig.hmftools.common.virus.AnnotatedVirus input,
            @NotNull VirusInterpreterEntry converted)
    {
        assertEquals(input.name(), converted.name());
        assertEquals(input.qcStatus().name(), converted.qcStatus().name());
        assertEquals(input.integrations(), converted.integrations());
        assertEquals(input.percentageCovered(), converted.percentageCovered(), EPSILON);
        assertEquals(input.meanCoverage(), converted.meanCoverage(), EPSILON);
        assertEquals(input.expectedClonalCoverage(), converted.expectedClonalCoverage());
        assertEquals(input.reported(), converted.reported());
        assertEquals(input.virusDriverLikelihoodType().name(), converted.driverLikelihood().name());
    }
}