package com.hartwig.hmftools.common.virus;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class VirusTypeTest
{
    @Test
    public void canExtractVirusConstants()
    {
        VirusType virusTypeHPV = VirusType.fromVirusName("HPV");
        assertEquals(VirusType.HPV, virusTypeHPV);

        VirusType virusTypeMCV = VirusType.fromVirusName("MCV");
        assertEquals(VirusType.MCV, virusTypeMCV);

        VirusType virusTypeEBV = VirusType.fromVirusName("EBV");
        assertEquals(VirusType.EBV, virusTypeEBV);

        VirusType virusTypeHBV = VirusType.fromVirusName("HBV");
        assertEquals(VirusType.HBV, virusTypeHBV);

        VirusType virusTypeHHV8 = VirusType.fromVirusName("HHV-8");
        assertEquals(VirusType.HHV8, virusTypeHHV8);
    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownVirusConstants()
    {
        //noinspection ResultOfMethodCallIgnored
        VirusType.fromVirusName("ABC");
    }
}