package com.hartwig.hmftools.orange;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_DESC;
import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import com.hartwig.hmftools.common.pipeline.PipelineToolDirectories;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionFile;
import com.hartwig.hmftools.common.rna.GeneExpressionFile;
import com.hartwig.hmftools.common.rna.GeneFusionFile;
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

    String ISOFOX_GENE_DISTRIBUTION_CSV = "isofox_gene_distribution";
    String ISOFOX_ALT_SJ_COHORT_CSV = "isofox_alt_sj_cohort";

    static void registerConfig(@NotNull ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(RNA_SAMPLE_ID, false, "(Optional) The RNA sample of the tumor sample for which ORANGE will run");

        configBuilder.addPath(ISOFOX_GENE_DISTRIBUTION_CSV, false, "(Optional) Path to isofox gene distribution CSV");
        configBuilder.addPath(ISOFOX_ALT_SJ_COHORT_CSV, false, "(Optional) Path to isofox alt SJ cohort CSV.");

        configBuilder.addPath(ISOFOX_DIR_CFG, false, ISOFOX_DIR_DESC);
    }

    @NotNull
    String rnaSampleId();

    @NotNull
    String isofoxGeneDistributionCsv();

    @NotNull
    String isofoxAltSjCohortCsv();

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
        boolean hasGeneDistribution = configBuilder.hasValue(ISOFOX_GENE_DISTRIBUTION_CSV);
        boolean hasAltSjCohortFreq = configBuilder.hasValue(ISOFOX_ALT_SJ_COHORT_CSV);

        boolean anyConfigPresent = hasRnaSampleId || hasIsofoxDir || hasGeneDistribution || hasAltSjCohortFreq;

        if(!anyConfigPresent)
        {
            LOGGER.info("RNA config not present, will continue without RNA configuration.");
            return null;
        }

        boolean allConfigPresent = hasRnaSampleId && hasIsofoxDir && hasGeneDistribution && hasAltSjCohortFreq;

        if(!allConfigPresent)
        {
            throw new IllegalArgumentException(String.format(
                    "RNA missing required config items: rnaSampleId(%s) isofoxDir(%s) geneCohort(%s) altSjCohort(%s)",
                    hasRnaSampleId,
                    hasIsofoxDir,
                    hasGeneDistribution,
                    hasAltSjCohortFreq));
        }

        String rnaSampleId = configBuilder.getValue(RNA_SAMPLE_ID);
        LOGGER.debug("RNA sample configured as {}", rnaSampleId);


        String geneDistributionFile = configBuilder.getValue(ISOFOX_GENE_DISTRIBUTION_CSV);
        String altSpliceJuncCohortFile = configBuilder.getValue(ISOFOX_ALT_SJ_COHORT_CSV);

        String tumorSampleId = configBuilder.getValue(OrangeConfig.TUMOR_SAMPLE_ID);

        String isofoxDir = pathResolver.resolveMandatoryToolDirectory(ISOFOX_DIR_CFG, defaultToolDirectories.isofoxDir());
        String geneDataFile = PathUtil.mandatoryPath(GeneExpressionFile.generateFilename(isofoxDir, tumorSampleId));
        String statisticsFile = PathUtil.mandatoryPath(RnaStatistics.generateFilename(isofoxDir, tumorSampleId));
        String altSpliceJuncFile = PathUtil.mandatoryPath(AltSpliceJunctionFile.generateFilename(isofoxDir, tumorSampleId));
        String fusionsFile = PathUtil.mandatoryPath(GeneFusionFile.generateFilename(isofoxDir, tumorSampleId));

        return ImmutableOrangeRnaConfig.builder()
                .rnaSampleId(rnaSampleId)
                .isofoxGeneDistributionCsv(geneDistributionFile)
                .isofoxAltSjCohortCsv(altSpliceJuncCohortFile)
                .isofoxSummaryCsv(statisticsFile)
                .isofoxGeneDataCsv(geneDataFile)
                .isofoxFusionCsv(fusionsFile)
                .isofoxAltSpliceJunctionCsv(altSpliceJuncFile)
                .build();
    }
}
