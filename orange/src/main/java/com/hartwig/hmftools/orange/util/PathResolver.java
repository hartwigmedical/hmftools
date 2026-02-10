package com.hartwig.hmftools.orange.util;

import static java.lang.String.format;

import static com.hartwig.hmftools.orange.util.PathUtil.mandatoryPath;

import java.io.File;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PathResolver
{
    private final ConfigBuilder mConfigBuilder;
    @Nullable
    private final String mPipelineSampleRootDir;
    @Nullable
    private final String mSampleDataDir;

    public PathResolver(final ConfigBuilder configBuilder, @Nullable final String pipelineSampleRootDir,
            @Nullable final String sampleDataDir)
    {
        mConfigBuilder = configBuilder;
        mPipelineSampleRootDir = pipelineSampleRootDir;
        mSampleDataDir = sampleDataDir;
    }

    public String resolveMandatoryToolDirectory(final String toolDirConfigKey, final String defaultPipelineToolDir)
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
    public String resolveOptionalToolDirectory(final String toolDirConfigKey, final String defaultPipelineToolDir)
    {
        if(mConfigBuilder.hasValue(toolDirConfigKey))
        {
            return mConfigBuilder.getValue(toolDirConfigKey);
        }

        if(mPipelineSampleRootDir != null)
        {
            return mPipelineSampleRootDir + File.separator + defaultPipelineToolDir;
        }

        return mSampleDataDir;
    }

    public String resolveMandatoryToolPlotsDirectory(final String toolPlotDirConfigKey, final String defaultPipelineToolDir)
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
    public String resolveOptionalToolPlotsDirectory(final String toolPlotDirConfigKey, final String defaultPipelineToolDir)
    {
        if(mConfigBuilder.hasValue(toolPlotDirConfigKey))
        {
            return mConfigBuilder.getValue(toolPlotDirConfigKey);
        }

        if(mPipelineSampleRootDir != null)
        {
            return mPipelineSampleRootDir + File.separator + defaultPipelineToolDir + File.separator + "plot";
        }

        return null;
    }
}
