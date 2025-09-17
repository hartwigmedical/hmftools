package com.hartwig.hmftools.orange.util;

import static java.lang.String.format;

import static com.hartwig.hmftools.orange.util.PathUtil.mandatoryPath;

import java.io.File;

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
}
