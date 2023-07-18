package com.hartwig.hmftools.healthchecker.runners;

import java.util.Optional;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class HealthCheckSampleConfiguration
{
    public final String SampleName;
    public final String WgsMetricsFile;
    public final String FlagstatFile;

    public HealthCheckSampleConfiguration(final String sampleName, final String wgsMetricsFile, final String flagstatFile)
    {
        SampleName = sampleName;
        WgsMetricsFile = wgsMetricsFile;
        FlagstatFile = flagstatFile;
    }

    public HealthCheckSampleConfiguration of(@Nullable String sampleName, @Nullable String wgsMetricsFile, @Nullable String flagstatFile)
    {
        if(sampleName != null)
        {
            return new HealthCheckSampleConfiguration(sampleName, notnullOrThrow(wgsMetricsFile), notnullOrThrow(flagstatFile));
        }

        return null;
    }

    @NotNull
    static String notnullOrThrow(@Nullable final String wgsMetricsFile)
    {
        return Optional.ofNullable(wgsMetricsFile).orElseThrow();
    }
}
