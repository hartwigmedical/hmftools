package com.hartwig.hmftools.orange;

import static com.hartwig.hmftools.common.metrics.GeneDepthFile.generateGeneCoverageFilename;
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

import java.nio.file.Files;
import java.nio.file.Paths;

import com.hartwig.hmftools.common.chord.ChordDataFile;
import com.hartwig.hmftools.common.cuppa.CuppaPredictions;
import com.hartwig.hmftools.common.metrics.BamFlagStats;
import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.common.peach.PeachGenotypeFile;
import com.hartwig.hmftools.common.pipeline.PipelineToolDirectories;
import com.hartwig.hmftools.common.redux.BqrFile;
import com.hartwig.hmftools.common.sigs.SignatureAllocationFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.virus.AnnotatedVirusFile;
import com.hartwig.hmftools.orange.util.PathResolver;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface OrangeRefConfig
{
    String REFERENCE_SAMPLE_ID = "reference_sample_id";

    String REFERENCE_REDUX_DIR_CFG = "ref_redux_dir";
    String REFERENCE_REDUX_DIR_DESC = "Path to Redux reference files";

    static void registerConfig(@NotNull ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(REFERENCE_SAMPLE_ID,
                false,
                "(Optional) The reference sample of the tumor sample for which ORANGE will run.");
        configBuilder.addPath(REFERENCE_REDUX_DIR_CFG, false, REFERENCE_REDUX_DIR_DESC);
        configBuilder.addPath(REF_METRICS_DIR_CFG, false, REF_METRICS_DIR_DESC);
        configBuilder.addPath(LINX_GERMLINE_DIR_CFG, false, LINX_GERMLINE_DIR_DESC);
        configBuilder.addPath(SAGE_GERMLINE_DIR_CFG, false, SAGE_GERMLINE_DIR_DESC);
        configBuilder.addPath(VIRUS_DIR_CFG, false, VIRUS_DIR_DESC);
        configBuilder.addPath(CHORD_DIR_CFG, false, CHORD_DIR_DESC);
        configBuilder.addPath(CUPPA_DIR_CFG, false, CUPPA_DIR_DESC);
        configBuilder.addPath(PEACH_DIR_CFG, false, PEACH_DIR_DESC);
        configBuilder.addPath(SIGS_DIR_CFG, false, SIGS_DIR_DESC);

    }

    @NotNull
    static OrangeRefConfig createConfig(
            final ConfigBuilder configBuilder, final PathResolver pathResolver, final PipelineToolDirectories defaultToolDirectories)
    {
        ImmutableOrangeRefConfig.Builder builder = ImmutableOrangeRefConfig.builder();
        String tumorSampleId = configBuilder.getValue(TUMOR_SAMPLE_ID);

        // Params required for WGS, Tumor only
        String virusDir = pathResolver.resolveOptionalToolDirectory(VIRUS_DIR_CFG, defaultToolDirectories.virusInterpreterDir());
        if(virusDir != null)
        {
            builder.annotatedVirusTsv(AnnotatedVirusFile.generateFileName(virusDir, tumorSampleId));
        }

        String chordDir = pathResolver.resolveMandatoryToolDirectory(CHORD_DIR_CFG, defaultToolDirectories.chordDir());
        builder.chordPredictionTxt(mandatoryPath(ChordDataFile.generateFilename(chordDir, tumorSampleId)));

        String sigsDir = pathResolver.resolveMandatoryToolDirectory(SIGS_DIR_CFG, defaultToolDirectories.sigsDir());
        builder.sigsAllocationTsv(mandatoryPath(SignatureAllocationFile.generateFilename(sigsDir, tumorSampleId)));

        // optionally required for WGS, adding Reference
        String refSampleId = configBuilder.getValue(REFERENCE_SAMPLE_ID);
        if(refSampleId != null)
        {
            LOGGER.debug("Ref sample has been configured as {}.", refSampleId);
            builder.referenceSampleId(refSampleId);

            if(configBuilder.hasValue(REFERENCE_REDUX_DIR_CFG))
            {
                String reduxDir = configBuilder.getValue(REFERENCE_REDUX_DIR_CFG);
                builder.refSampleBqrPlot(BqrFile.generatePlotFilename(reduxDir, refSampleId));
            }

            String refMetricsDir = pathResolver.resolveMandatoryToolDirectory(
                    REF_METRICS_DIR_CFG, defaultToolDirectories.germlineMetricsDir());
            String geneCoverageFile = generateGeneCoverageFilename(refMetricsDir, refSampleId);

            if(Files.exists(Paths.get(geneCoverageFile)))
            {
                builder.germlineGeneCoverageTsv(geneCoverageFile);

            }
            else
            {
                String sageGermlineDir = pathResolver.resolveMandatoryToolDirectory(SAGE_GERMLINE_DIR_CFG, defaultToolDirectories.sageGermlineDir());
                String legacySageGermlineCoverageFile = generateGeneCoverageFilenameLegacySage(sageGermlineDir, refSampleId);
                builder.germlineGeneCoverageTsv(legacySageGermlineCoverageFile);
            }

            String linxGermlineDir = pathResolver.resolveMandatoryToolDirectory(LINX_GERMLINE_DIR_CFG, defaultToolDirectories.linxGermlineDir());
            builder.linxGermlineDataDirectory(linxGermlineDir);
            
            String cuppaDir = pathResolver.resolveOptionalToolDirectory(CUPPA_DIR_CFG, defaultToolDirectories.cuppaDir());
            if(cuppaDir != null)
            {
                builder.cuppaVisDataTsv(mandatoryPath(CuppaPredictions.generateVisDataTsvFilename(cuppaDir, tumorSampleId)));
                builder.cuppaSummaryPlot(mandatoryPath(CuppaPredictions.generateVisPlotFilename(cuppaDir, tumorSampleId)));
            }

            // PEACH optional so that skipping it in oncoanalyser still generates an ORANGE report
            String peachDir = pathResolver.resolveOptionalToolDirectory(PEACH_DIR_CFG, defaultToolDirectories.peachDir());
            if(peachDir != null)
            {
                String peachGenotypeTsv = mandatoryPath(PeachGenotypeFile.generateFileName(peachDir, refSampleId));
                builder.peachGenotypeTsv(peachGenotypeTsv);
            }

            builder.refSampleWGSMetricsFile(mandatoryPath(BamMetricSummary.generateFilename(refMetricsDir, refSampleId)));
            builder.refSampleFlagstatFile(mandatoryPath(BamFlagStats.generateFilename(refMetricsDir, refSampleId)));
        }

        return builder.build();
    }

    private static String generateGeneCoverageFilenameLegacySage(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + ".sage.gene.coverage.tsv";
    }

    @Nullable
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
    String germlineGeneCoverageTsv();

    @Nullable
    String refSampleBqrPlot();

    @Nullable
    String linxGermlineDataDirectory();

    @Nullable
    String peachGenotypeTsv();

    @Nullable
    String refSampleWGSMetricsFile();

    @Nullable
    String refSampleFlagstatFile();
}
