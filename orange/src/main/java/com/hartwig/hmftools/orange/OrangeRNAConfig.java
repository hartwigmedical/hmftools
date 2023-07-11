package com.hartwig.hmftools.orange;

import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.ISOFOX_DIR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PIPELINE_SAMPLE_ROOT_DIR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.orange.OrangeConfig.getToolDirectory;

import com.hartwig.hmftools.common.rna.AltSpliceJunctionFile;
import com.hartwig.hmftools.common.rna.GeneExpressionFile;
import com.hartwig.hmftools.common.rna.GeneFusionFile;
import com.hartwig.hmftools.common.rna.RnaStatistics;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.orange.util.Config;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface OrangeRNAConfig {

    Logger LOGGER = LogManager.getLogger(OrangeRNAConfig.class);

    String RNA_SAMPLE_ID = "rna_sample_id";

    String ISOFOX_GENE_DISTRIBUTION_CSV = "isofox_gene_distribution";
    String ISOFOX_ALT_SJ_COHORT_CSV = "isofox_alt_sj_cohort";

    static void registerConfig(final ConfigBuilder configBuilder) {

        configBuilder.addConfigItem(RNA_SAMPLE_ID, "(Optional) The RNA sample of the tumor sample for which ORANGE will run");

        configBuilder.addPath(ISOFOX_GENE_DISTRIBUTION_CSV, false, "(Optional) Path to isofox gene distribution CSV");
        configBuilder.addPath(ISOFOX_ALT_SJ_COHORT_CSV, false, "(Optional) Path to isofox alt SJ cohort CSV.");

        configBuilder.addPath(ISOFOX_DIR_CFG, false, ISOFOX_DIR_DESC);
    }

    String rnaSampleId();

    String isofoxGeneDistributionCsv();

    String isofoxAltSjCohortCsv();

    String isofoxSummaryCsv();

    String isofoxGeneDataCsv();

    String isofoxFusionCsv();

    String isofoxAltSpliceJunctionCsv();

    static OrangeRNAConfig createConfig(final ConfigBuilder configBuilder) {

        boolean hasRnaSampleId = configBuilder.hasValue(RNA_SAMPLE_ID);
        boolean hasIsofoxDir = configBuilder.hasValue(ISOFOX_DIR_CFG);
        boolean hasGeneDistribution = configBuilder.hasValue(ISOFOX_GENE_DISTRIBUTION_CSV);
        boolean hasAltSjCohortFreq = configBuilder.hasValue(ISOFOX_ALT_SJ_COHORT_CSV);

        boolean anyConfigPresent = hasRnaSampleId || hasIsofoxDir || hasGeneDistribution || hasAltSjCohortFreq;

        if (!anyConfigPresent) {
            LOGGER.info("RNA config not present, will continue without RNA configuration.");
            return null;
        }

        boolean allConfigPresent = hasRnaSampleId && hasIsofoxDir && hasGeneDistribution && hasAltSjCohortFreq;

        if (!allConfigPresent) {
            throw new IllegalArgumentException(String.format(
                    "RNA missing required config items: rnaSampleId(%s) isofoxDir(%s) geneCohort(%s) altSjCohort(%s)",
                    hasRnaSampleId,
                    hasIsofoxDir,
                    hasGeneDistribution,
                    hasAltSjCohortFreq));
        }

        String rnaSampleId = configBuilder.getValue(RNA_SAMPLE_ID);
        LOGGER.debug("RNA sample configured as {}", rnaSampleId);

        String pipelineSampleRootDir = checkAddDirSeparator(configBuilder.getValue(PIPELINE_SAMPLE_ROOT_DIR));
        String sampleDataDir = checkAddDirSeparator(configBuilder.getValue(SAMPLE_DATA_DIR_CFG));

        String isofoxDir = getToolDirectory(configBuilder, pipelineSampleRootDir, sampleDataDir, ISOFOX_DIR_CFG, ISOFOX_DIR);

        String tumorSampleId = configBuilder.getValue(OrangeConfig.TUMOR_SAMPLE_ID);
        String geneDataFile = Config.fileExists(GeneExpressionFile.generateFilename(isofoxDir, tumorSampleId));
        String statisticsFile = Config.fileExists(RnaStatistics.generateFilename(isofoxDir, tumorSampleId));
        String altSpliceJuncFile = Config.fileExists(AltSpliceJunctionFile.generateFilename(isofoxDir, tumorSampleId));
        String fusionsFile = Config.fileExists(GeneFusionFile.generateFilename(isofoxDir, tumorSampleId));
        String geneDistributionFile = Config.fileExists(configBuilder.getValue(ISOFOX_GENE_DISTRIBUTION_CSV));
        String altSpliceJuncCohortFile = Config.fileExists(configBuilder.getValue(ISOFOX_ALT_SJ_COHORT_CSV));

        return ImmutableOrangeRNAConfig.builder()
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
