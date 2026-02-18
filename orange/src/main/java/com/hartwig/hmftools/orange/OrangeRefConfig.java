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
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_DESC;
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

public class OrangeRefConfig
{
    public final String ReferenceId;

    public final String AnnotatedVirusTsv;
    public final String ChordPredictionTxt;
    public final String CuppaVisDataTsv;
    public final String CuppaSummaryPlot;
    public final String SigsAllocationTsv;
    public final String GermlineGeneCoverageTsv;
    public final String RefSampleBqrPlot;
    public final String LinxGermlineDataDirectory;
    public final String PeachGenotypeTsv;
    public final String RefSampleWGSMetricsFile;
    public final String RefSampleFlagstatFile;

    private static String REFERENCE_REDUX_DIR_CFG = "ref_redux_dir";
    private static String REFERENCE_REDUX_DIR_DESC = "Path to Redux reference files";

    public OrangeRefConfig(
            final ConfigBuilder configBuilder, final String tumorSampleId, final PathResolver pathResolver,
            final PipelineToolDirectories defaultToolDirectories)
    {
        ReferenceId = configBuilder.getValue(REFERENCE);

        if(ReferenceId != null)
        {
            LOGGER.debug("Ref sample has been configured as {}.", ReferenceId);

            String reduxDir = configBuilder.getValue(REFERENCE_REDUX_DIR_CFG);
            RefSampleBqrPlot = BqrFile.generatePlotFilename(reduxDir, ReferenceId);

            String refMetricsDir = pathResolver.resolveMandatoryToolDirectory(
                    REF_METRICS_DIR_CFG, defaultToolDirectories.germlineMetricsDir());
            String geneCoverageFile = generateGeneCoverageFilename(refMetricsDir, ReferenceId);

            if(Files.exists(Paths.get(geneCoverageFile)))
            {
                GermlineGeneCoverageTsv = geneCoverageFile;
            }
            else
            {
                GermlineGeneCoverageTsv = "";
            }

            String linxGermlineDir = pathResolver.resolveMandatoryToolDirectory(LINX_GERMLINE_DIR_CFG, defaultToolDirectories.linxGermlineDir());
            LinxGermlineDataDirectory = linxGermlineDir;

            String cuppaDir = pathResolver.resolveOptionalToolDirectory(CUPPA_DIR_CFG, defaultToolDirectories.cuppaDir());
            if(cuppaDir != null)
            {
                CuppaVisDataTsv = mandatoryPath(CuppaPredictions.generateVisDataTsvFilename(cuppaDir, tumorSampleId));
                CuppaSummaryPlot = mandatoryPath(CuppaPredictions.generateVisPlotFilename(cuppaDir, tumorSampleId));
            }
            else
            {
                CuppaVisDataTsv = "";
                CuppaSummaryPlot = "";
            }

            // PEACH optional so that skipping it in oncoanalyser still generates an ORANGE report
            String peachDir = pathResolver.resolveOptionalToolDirectory(PEACH_DIR_CFG, defaultToolDirectories.peachDir());
            if(peachDir != null)
            {
                String peachGenotypeTsv = mandatoryPath(PeachGenotypeFile.generateFileName(peachDir, ReferenceId));
                PeachGenotypeTsv = peachGenotypeTsv;
            }
            else
            {
                PeachGenotypeTsv = "";
            }

            RefSampleWGSMetricsFile = mandatoryPath(BamMetricSummary.generateFilename(refMetricsDir, ReferenceId));
            RefSampleFlagstatFile = mandatoryPath(BamFlagStats.generateFilename(refMetricsDir, ReferenceId));

            String chordDir = pathResolver.resolveMandatoryToolDirectory(CHORD_DIR_CFG, defaultToolDirectories.chordDir());
            ChordPredictionTxt = mandatoryPath(ChordDataFile.generateFilename(chordDir, tumorSampleId));

            String sigsDir = pathResolver.resolveMandatoryToolDirectory(SIGS_DIR_CFG, defaultToolDirectories.sigsDir());
            SigsAllocationTsv = mandatoryPath(SignatureAllocationFile.generateFilename(sigsDir, tumorSampleId));

            String virusDir = pathResolver.resolveOptionalToolDirectory(VIRUS_DIR_CFG, defaultToolDirectories.virusInterpreterDir());
            if(virusDir != null)
            {
                AnnotatedVirusTsv = AnnotatedVirusFile.generateFileName(virusDir, tumorSampleId);
            }
            else
            {
                AnnotatedVirusTsv = "";
            }
        }
        else
        {
            AnnotatedVirusTsv = "";
            ChordPredictionTxt = "";
            CuppaVisDataTsv = "";
            CuppaSummaryPlot = "";
            SigsAllocationTsv = "";
            GermlineGeneCoverageTsv = "";
            RefSampleBqrPlot = "";
            LinxGermlineDataDirectory = "";
            PeachGenotypeTsv = "";
            RefSampleWGSMetricsFile = "";
            RefSampleFlagstatFile = "";
        }
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(REFERENCE, false, REFERENCE_DESC);
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
}
