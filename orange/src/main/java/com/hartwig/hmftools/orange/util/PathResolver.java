package com.hartwig.hmftools.orange.util;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.METRICS_DIR;
import static com.hartwig.hmftools.orange.util.PathUtil.mandatoryPath;

import java.io.File;

import com.hartwig.hmftools.common.metrics.BamFlagStats;
import com.hartwig.hmftools.common.metrics.BamMetricsSummary;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PathResolver
{
    @NotNull
    private final ConfigBuilder configBuilder;
    @Nullable
    private final String pipelineSampleRootDir;
    @Nullable
    private final String sampleDataDir;

    public PathResolver(@NotNull final ConfigBuilder configBuilder, @Nullable final String pipelineSampleRootDir,
            @Nullable final String sampleDataDir)
    {
        this.configBuilder = configBuilder;
        this.pipelineSampleRootDir = pipelineSampleRootDir;
        this.sampleDataDir = sampleDataDir;
    }

    @NotNull
    public String resolveMandatoryToolDirectory(@NotNull String toolDirConfigKey, @NotNull String defaultPipelineToolDir)
    {
        String toolDir = resolveOptionalToolDirectory(toolDirConfigKey, defaultPipelineToolDir);
        if(toolDir == null)
        {
            throw new IllegalArgumentException(format("Failed to determine tool directory for configuration [%s/%s].",
                    toolDirConfigKey, defaultPipelineToolDir));
        }
        return mandatoryPath(toolDir);
    }

    @Nullable
    public String resolveOptionalToolDirectory(@NotNull String toolDirConfigKey, @NotNull String defaultPipelineToolDir)
    {
        if(configBuilder.hasValue(toolDirConfigKey))
        {
            return configBuilder.getValue(toolDirConfigKey);
        }

        if(pipelineSampleRootDir != null)
        {
            return pipelineSampleRootDir + File.separator + defaultPipelineToolDir;
        }

        return sampleDataDir;
    }

    @NotNull
    public String resolveMandatoryToolPlotsDirectory(@NotNull String toolPlotDirConfigKey, @NotNull String defaultPipelineToolDir)
    {
        String plotDir = resolveOptionalToolPlotsDirectory(toolPlotDirConfigKey, defaultPipelineToolDir);
        if(plotDir == null)
        {
            throw new IllegalArgumentException(format("Failed to determine plot directory for configuration [%s/%s].",
                    toolPlotDirConfigKey, defaultPipelineToolDir));
        }
        return mandatoryPath(plotDir);
    }

    @Nullable
    public String resolveOptionalToolPlotsDirectory(@NotNull String toolPlotDirConfigKey, @NotNull String defaultPipelineToolDir)
    {
        if(configBuilder.hasValue(toolPlotDirConfigKey))
        {
            return configBuilder.getValue(toolPlotDirConfigKey);
        }

        if(pipelineSampleRootDir != null)
        {
            return pipelineSampleRootDir + File.separator + defaultPipelineToolDir + File.separator + "plot";
        }

        return null;
    }

    @NotNull
    public String resolveMandatoryMetricsFile(@NotNull String metricFileConfigKey, @NotNull String defaultPipelineToolDir,
            @NotNull String sampleId)
    {
        String metricsFile = resolveOptionalMetricsFile(metricFileConfigKey, defaultPipelineToolDir, sampleId);
        if(metricsFile == null)
        {
            throw new IllegalArgumentException(format("Failed to determine metrics file for configuration [%s/%s].",
                    metricFileConfigKey, defaultPipelineToolDir));
        }
        return mandatoryPath(metricsFile);
    }

    @Nullable
    public String resolveOptionalMetricsFile(@NotNull String metricFileConfigKey, @NotNull String defaultPipelineToolDir,
            @NotNull String sampleId)
    {
        if(configBuilder.hasValue(metricFileConfigKey))
        {
            return configBuilder.getValue(metricFileConfigKey);
        }

        String directory = resolveMetricsDirectory(defaultPipelineToolDir, sampleId);
        if(directory == null)
        {
            return null;
        }

        return defaultPipelineToolDir.equals(METRICS_DIR) ?
                BamMetricsSummary.generateFilename(directory, sampleId) : BamFlagStats.generateFilename(directory, sampleId);
    }

    @Nullable
    private String resolveMetricsDirectory(@NotNull String defaultPipelineToolDir, @NotNull String sampleId)
    {
        if(pipelineSampleRootDir == null)
        {
            return null;
        }

        return pipelineSampleRootDir + File.separator + sampleId + File.separator + defaultPipelineToolDir;
    }
}
