package com.hartwig.hmftools.orange;

import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.SAGE_GERMLINE_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.SAGE_SOMATIC_DIR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PIPELINE_SAMPLE_ROOT_DIR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_GERMLINE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.orange.OrangeConfig.REFERENCE_SAMPLE_ID;
import static com.hartwig.hmftools.orange.OrangeConfig.getToolDirectory;
import static com.hartwig.hmftools.orange.util.Config.fileIfExists;

import com.hartwig.hmftools.common.sage.SageCommon;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface OrangeRefConfig
{
    @NotNull
    String sageGermlineGeneCoverageTsv();

    @NotNull
    String sageSomaticRefSampleBQRPlot();

    @Nullable
    static OrangeRefConfig createConfig(@NotNull ConfigBuilder configBuilder, @NotNull String sageToolDirectory)
    {
        ImmutableOrangeRefConfig.Builder builder = ImmutableOrangeRefConfig.builder();
        String refSampleId = configBuilder.getValue(REFERENCE_SAMPLE_ID);

        if(refSampleId != null)
        {
            return null;
        }

        // TODO these config directories are also extracted elsewhere, maybe share?
        String pipelineSampleRootDir = checkAddDirSeparator(configBuilder.getValue(PIPELINE_SAMPLE_ROOT_DIR));
        String sampleDataDir = checkAddDirSeparator(configBuilder.getValue(SAMPLE_DATA_DIR_CFG));
        String sageSomaticDir = getToolDirectory(configBuilder, pipelineSampleRootDir, sampleDataDir, SAGE_DIR_CFG, SAGE_SOMATIC_DIR);

        // TODO sageGermline should be here or in WGS?
        String sageGermlineDir = getToolDirectory(
                configBuilder, pipelineSampleRootDir, sampleDataDir, SAGE_GERMLINE_DIR_CFG, SAGE_GERMLINE_DIR);

        builder.sageGermlineGeneCoverageTsv(fileIfExists(SageCommon.generateGeneCoverageFilename(sageGermlineDir, refSampleId)));

        // TODO why are somatic dir and ref sample id together here?
        builder.sageSomaticRefSampleBQRPlot(fileIfExists(SageCommon.generateBqrPlotFilename(sageSomaticDir, refSampleId)));

        return builder.build();
    }
}
