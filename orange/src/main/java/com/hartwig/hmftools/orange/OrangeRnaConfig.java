package com.hartwig.hmftools.orange;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_DESC;
import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import com.hartwig.hmftools.common.pipeline.PipelineToolDirectories;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionFile;
import com.hartwig.hmftools.common.rna.GeneExpressionFile;
import com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile;
import com.hartwig.hmftools.common.rna.RnaFusionFile;
import com.hartwig.hmftools.common.rna.RnaStatisticFile;
import com.hartwig.hmftools.common.rna.RnaStatistics;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.orange.util.PathResolver;
import com.hartwig.hmftools.orange.util.PathUtil;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface OrangeRnaConfig
{
    String RNA_SAMPLE_ID = "rna_sample_id";

    static void registerConfig(@NotNull ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(RNA_SAMPLE_ID, false, "(Optional) The RNA sample of the tumor sample for which ORANGE will run");
        configBuilder.addPath(ISOFOX_DIR_CFG, false, ISOFOX_DIR_DESC);
    }

    @NotNull
    String rnaSampleId();

    @NotNull
    String isofoxSummaryCsv();

    @NotNull
    String isofoxGeneDataCsv();

    @NotNull
    String isofoxFusionCsv();

    @NotNull
    String isofoxAltSpliceJunctionCsv();

    @Nullable
    static OrangeRnaConfig createConfig(
            final ConfigBuilder configBuilder, final PathResolver pathResolver, final PipelineToolDirectories defaultToolDirectories)
    {
        boolean hasRnaSampleId = configBuilder.hasValue(RNA_SAMPLE_ID);
        boolean hasIsofoxDir = configBuilder.hasValue(ISOFOX_DIR_CFG);

        boolean anyConfigPresent = hasRnaSampleId || hasIsofoxDir;

        if(!anyConfigPresent)
        {
            LOGGER.info("RNA config not present, will continue without RNA configuration.");
            return null;
        }

        boolean allConfigPresent = hasRnaSampleId && hasIsofoxDir;

        if(!allConfigPresent)
        {
            throw new IllegalArgumentException(String.format(
                    "RNA missing required config items: rnaSampleId(%s) isofoxDir(%s)",
                    hasRnaSampleId, hasIsofoxDir));
        }

        String rnaSampleId = configBuilder.getValue(RNA_SAMPLE_ID);
        LOGGER.debug("RNA sample configured as {}", rnaSampleId);

        String tumorSampleId = configBuilder.getValue(OrangeConfig.TUMOR_SAMPLE_ID);

        String isofoxDir = pathResolver.resolveMandatoryToolDirectory(ISOFOX_DIR_CFG, defaultToolDirectories.isofoxDir());
        String geneDataFile = PathUtil.mandatoryPath(GeneExpressionFile.generateFilename(isofoxDir, tumorSampleId));
        String statisticsFile = PathUtil.mandatoryPath(RnaStatisticFile.generateFilename(isofoxDir, tumorSampleId));
        String altSpliceJuncFile = PathUtil.mandatoryPath(NovelSpliceJunctionFile.generateFilename(isofoxDir, tumorSampleId));
        String fusionsFile = PathUtil.mandatoryPath(RnaFusionFile.generateFilename(isofoxDir, tumorSampleId));

        return ImmutableOrangeRnaConfig.builder()
                .rnaSampleId(rnaSampleId)
                .isofoxSummaryCsv(statisticsFile)
                .isofoxGeneDataCsv(geneDataFile)
                .isofoxFusionCsv(fusionsFile)
                .isofoxAltSpliceJunctionCsv(altSpliceJuncFile)
                .build();
    }
}
