package com.hartwig.hmftools.healthchecker.runners;

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
}
