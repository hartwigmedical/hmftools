package com.hartwig.hmftools.orange;

import static com.hartwig.hmftools.common.metrics.GeneDepthFile.generateGeneCoverageFilename;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.CHORD_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.CUPPA_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.LINX_GERMLINE_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.METRICS_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.PEACH_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.SAGE_GERMLINE_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.SAGE_SOMATIC_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.SIGS_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.VIRUS_INTERPRETER_DIR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CHORD_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CHORD_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CUPPA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CUPPA_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_GERMLINE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_GERMLINE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PEACH_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PEACH_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REF_METRICS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REF_METRICS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_GERMLINE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_GERMLINE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SIGS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SIGS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.VIRUS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.VIRUS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;
import static com.hartwig.hmftools.orange.OrangeConfig.TUMOR_SAMPLE_ID;
import static com.hartwig.hmftools.orange.util.PathUtil.mandatoryPath;

import com.hartwig.hmftools.common.chord.ChordDataFile;
import com.hartwig.hmftools.common.cuppa.CuppaPredictions;
import com.hartwig.hmftools.common.metrics.BamFlagStats;
import com.hartwig.hmftools.common.metrics.BamMetricsSummary;
import com.hartwig.hmftools.common.peach.PeachGenotypeFile;
import com.hartwig.hmftools.common.sage.SageCommon;
import com.hartwig.hmftools.common.sigs.SignatureAllocationFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.virus.AnnotatedVirusFile;
import com.hartwig.hmftools.orange.util.PathResolver;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface OrangeWGSRefConfig {
    String REFERENCE_SAMPLE_ID = "reference_sample_id";

    static void registerConfig(@NotNull ConfigBuilder configBuilder) {
        configBuilder.addConfigItem(REFERENCE_SAMPLE_ID,
                false,
                "(Optional) The reference sample of the tumor sample for which ORANGE will run.");
        configBuilder.addPath(REF_METRICS_DIR_CFG, false, REF_METRICS_DIR_DESC);
        configBuilder.addPath(LINX_GERMLINE_DIR_CFG, false, LINX_GERMLINE_DIR_DESC);
        configBuilder.addPath(SAGE_GERMLINE_DIR_CFG, false, SAGE_GERMLINE_DIR_DESC);
        configBuilder.addPath(VIRUS_DIR_CFG, false, VIRUS_DIR_DESC);
        configBuilder.addPath(CHORD_DIR_CFG, false, CHORD_DIR_DESC);
        configBuilder.addPath(CUPPA_DIR_CFG, false, CUPPA_DIR_DESC);
        configBuilder.addPath(PEACH_DIR_CFG, false, PEACH_DIR_DESC);
        configBuilder.addPath(SIGS_DIR_CFG, false, SIGS_DIR_DESC);

    }

    // params for WGS Tumor only
    @NotNull
    String annotatedVirusTsv();

    @NotNull
    String chordPredictionTxt();

    @Nullable
    String cuppaVisDataTsv();

    @Nullable
    String cuppaSummaryPlot();

    @NotNull
    String sigsAllocationTsv();

    // additional params for WGS Ref
    @Nullable
    String referenceSampleId();

    @Nullable
    String sageGermlineGeneCoverageTsv();

    @Nullable
    String sageSomaticRefSampleBQRPlot();

    @Nullable
    String linxGermlineDataDirectory();

    @Nullable
    String peachGenotypeTsv();

    @Nullable
    String refSampleWGSMetricsFile();

    @Nullable
    String refSampleFlagstatFile();

    @NotNull
    static OrangeWGSRefConfig createConfig(@NotNull ConfigBuilder configBuilder, @NotNull PathResolver pathResolver) {
        ImmutableOrangeWGSRefConfig.Builder builder = ImmutableOrangeWGSRefConfig.builder();
        String tumorSampleId = configBuilder.getValue(TUMOR_SAMPLE_ID);

        // Params required for WGS, Tumor only
        String virusDir = pathResolver.resolveMandatoryToolDirectory(VIRUS_DIR_CFG, VIRUS_INTERPRETER_DIR);
        builder.annotatedVirusTsv(mandatoryPath(AnnotatedVirusFile.generateFileName(virusDir, tumorSampleId)));

        String chordDir = pathResolver.resolveMandatoryToolDirectory(CHORD_DIR_CFG, CHORD_DIR);
        builder.chordPredictionTxt(mandatoryPath(ChordDataFile.generateFilename(chordDir, tumorSampleId)));

        String sigsDir = pathResolver.resolveMandatoryToolDirectory(SIGS_DIR_CFG, SIGS_DIR);
        builder.sigsAllocationTsv(mandatoryPath(SignatureAllocationFile.generateFilename(sigsDir, tumorSampleId)));

        // optionally required for WGS, adding Reference
        String refSampleId = configBuilder.getValue(REFERENCE_SAMPLE_ID);
        if (refSampleId != null) {
            LOGGER.debug("Ref sample has been configured as {}.", refSampleId);
            builder.referenceSampleId(refSampleId);

            String sageSomaticDir = pathResolver.resolveMandatoryToolDirectory(SAGE_DIR_CFG, SAGE_SOMATIC_DIR);
            builder.sageSomaticRefSampleBQRPlot(mandatoryPath(SageCommon.generateBqrPlotFilename(sageSomaticDir, refSampleId)));

            String sageGermlineDir = pathResolver.resolveMandatoryToolDirectory(SAGE_GERMLINE_DIR_CFG, SAGE_GERMLINE_DIR);
            builder.sageGermlineGeneCoverageTsv(backwardsCompatibleSageGermlineGeneConverage(sageGermlineDir, refSampleId));

            String linxGermlineDir = pathResolver.resolveMandatoryToolDirectory(LINX_GERMLINE_DIR_CFG, LINX_GERMLINE_DIR);
            builder.linxGermlineDataDirectory(linxGermlineDir);

            String cuppaDir = pathResolver.resolveMandatoryToolDirectory(CUPPA_DIR_CFG, CUPPA_DIR);
            builder.cuppaVisDataTsv(mandatoryPath(CuppaPredictions.generateVisDataTsvFilename(cuppaDir, tumorSampleId)));
            builder.cuppaSummaryPlot(mandatoryPath(CuppaPredictions.generateVisPlotFilename(cuppaDir, tumorSampleId)));

            // PEACH optional so that skipping it in oncoanalyser still generates an ORANGE report
            String peachDir = pathResolver.resolveOptionalToolDirectory(PEACH_DIR_CFG, PEACH_DIR);
            if (peachDir != null) {
                String peachGenotypeTsv = mandatoryPath(PeachGenotypeFile.generateFileName(peachDir, refSampleId));
                builder.peachGenotypeTsv(peachGenotypeTsv);
            }

            String refMetricsDir = pathResolver.resolveMandatoryToolDirectory(REF_METRICS_DIR_CFG, METRICS_DIR);
            builder.refSampleWGSMetricsFile(mandatoryPath(BamMetricsSummary.generateFilename(refMetricsDir, refSampleId)));
            builder.refSampleFlagstatFile(mandatoryPath(BamFlagStats.generateFilename(refMetricsDir, refSampleId)));
        }

        return builder.build();
    }

    @NotNull
    private static String backwardsCompatibleSageGermlineGeneConverage(String sageGermlineDir, String refSampleId) {
        try {
            String sageLegacyPath = generateGeneCoverageFilenameLegacySage(sageGermlineDir, refSampleId);
            LOGGER.info("Locating sage germline gene coverage [{}]", sageLegacyPath);
            return mandatoryPath(sageLegacyPath);
        } catch (IllegalArgumentException e) {
            return mandatoryPath(generateGeneCoverageFilename(sageGermlineDir, refSampleId));
        }
    }

    private static String generateGeneCoverageFilenameLegacySage(final String basePath, final String sample) {
        return checkAddDirSeparator(basePath) + sample + ".sage.gene.coverage.tsv";
    }
}
