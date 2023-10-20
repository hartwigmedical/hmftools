package com.hartwig.hmftools.orange.util;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.METRICS_DIR;
import static com.hartwig.hmftools.orange.util.Config.fileIfExists;

import java.io.File;

import com.hartwig.hmftools.common.flagstat.FlagstatFile;
import com.hartwig.hmftools.common.metrics.WGSMetricsFile;
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

    @Nullable
    public String resolveToolDirectory(@NotNull String toolDirConfigKey, @NotNull String defaultPipelineToolDir)
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
    public String resolveToolDirectoryIfExists(@NotNull String toolDirConfigKey, @NotNull String defaultPipelineToolDir)
    {
        String toolDir = resolveToolDirectory(toolDirConfigKey, defaultPipelineToolDir);
        if(toolDir == null)
        {
            throw new IllegalArgumentException(format("Failed to determine tool directory for configuration [%s/%s].", toolDirConfigKey, defaultPipelineToolDir));
        }
        return fileIfExists(toolDir);
    }

    @NotNull
    public String resolveToolPlotsDirectory(@NotNull String toolPlotDirConfigKey, @NotNull String defaultPipelineToolDir)
    {
        if(configBuilder.hasValue(toolPlotDirConfigKey))
        {
            return configBuilder.getValue(toolPlotDirConfigKey);
        }

        String plotDir = pipelineSampleRootDir != null ? fileIfExists(
                pipelineSampleRootDir + File.separator + defaultPipelineToolDir + File.separator + "plot") : null;
        if(plotDir == null)
        {
            throw new IllegalArgumentException(
                    "Plot directory cannot be determined. Please define either the tool directory or the sample directory.");
        }
        return plotDir;
    }

    @NotNull
    public String resolveMetricsFile(@NotNull String metricFileConfigKey, @NotNull String defaultPipelineToolDir, @NotNull String sampleId)
    {
        if(configBuilder.hasValue(metricFileConfigKey))
        {
            return configBuilder.getValue(metricFileConfigKey);
        }

        String directory = resolveMetricsDirectory(defaultPipelineToolDir, sampleId);
        if(directory == null)
        {
            throw new IllegalArgumentException(
                    "Metrics directory cannot be determined. Please define either the tool directory or the sample directory.");
        }

        return defaultPipelineToolDir.equals(METRICS_DIR) ?
                WGSMetricsFile.generateFilename(directory, sampleId) : FlagstatFile.generateFilename(directory, sampleId);
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
