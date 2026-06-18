package com.hartwig.hmftools.orange;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.time.LocalDate;

import org.junit.Test;

public class OrangeConfigTest
{
    @Test
    public void samplingDateIsNullWhenNotProvided()
    {
        OrangeConfig config = TestOrangeConfigFactory.createMinimalConfig();

        assertNull(config.SamplingDate);
    }

    @Test
    public void samplingDateParsesValidDate()
    {
        LocalDate result = OrangeConfig.interpretSamplingDateParam("210101");

        assertEquals(LocalDate.of(2021, 1, 1), result);
    }

    @Test
    public void samplingDateIsNullOnInvalidInput()
    {
        LocalDate result = OrangeConfig.interpretSamplingDateParam("not-a-date");

        assertNull(result);
    }

    @Test
    public void analysisDateDefaultsToTodayWhenNotProvided()
    {
        OrangeConfig config = TestOrangeConfigFactory.createMinimalConfig();

        assertEquals(LocalDate.now(), config.AnalysisDate);
    }

    @Test
    public void analysisDateParsesValidDate()
    {
        LocalDate result = OrangeConfig.interpretAnalysisDateParam("210101");

        assertEquals(LocalDate.of(2021, 1, 1), result);
    }

    @Test
    public void analysisDateDefaultsToTodayOnInvalidInput()
    {
        LocalDate result = OrangeConfig.interpretAnalysisDateParam("not-a-date");

        assertEquals(LocalDate.now(), result);
    }
}