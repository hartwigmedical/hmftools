package com.hartwig.hmftools.orange;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.time.LocalDate;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

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
        ConfigBuilder configBuilder = new ConfigBuilder();
        configBuilder.addConfigItem("sampling_date", false, "");
        configBuilder.parseCommandLine(new String[] { "-sampling_date", "210101" });

        assertEquals(LocalDate.of(2021, 1, 1), OrangeConfig.interpretDateParam(configBuilder, "sampling_date", null));
    }

    @Test
    public void samplingDateIsNullOnInvalidInput()
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        configBuilder.addConfigItem("sampling_date", false, "");
        configBuilder.parseCommandLine(new String[] { "-sampling_date", "not-a-date" });

        assertNull(OrangeConfig.interpretDateParam(configBuilder, "sampling_date", null));
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
        ConfigBuilder configBuilder = new ConfigBuilder();
        configBuilder.addConfigItem("analysis_date", false, "");
        configBuilder.parseCommandLine(new String[] { "-analysis_date", "210101" });

        assertEquals(LocalDate.of(2021, 1, 1), OrangeConfig.interpretDateParam(configBuilder, "analysis_date", LocalDate.now()));
    }

    @Test
    public void analysisDateDefaultsToTodayOnInvalidInput()
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        configBuilder.addConfigItem("analysis_date", false, "");
        configBuilder.parseCommandLine(new String[] { "-analysis_date", "not-a-date" });

        assertEquals(LocalDate.now(), OrangeConfig.interpretDateParam(configBuilder, "analysis_date", LocalDate.now()));
    }
}
