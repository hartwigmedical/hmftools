package com.hartwig.hmftools.isofox.cohort;

import static com.hartwig.hmftools.common.rna.RnaStatisticFile.SUMMARY_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.APP_NAME;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.SUMMARY;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.results.SummaryStats.loadFile;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.rna.RnaStatisticFile;
import com.hartwig.hmftools.common.rna.RnaStatistics;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.isofox.expression.cohort.ExpressionCohortCompare;
import com.hartwig.hmftools.isofox.expression.cohort.ExpressionCohortDistribution;
import com.hartwig.hmftools.isofox.expression.cohort.ExpressionMatrix;
import com.hartwig.hmftools.isofox.expression.cohort.ExternalExpressionCompare;
import com.hartwig.hmftools.isofox.expression.cohort.GeneratePanelNormalisation;
import com.hartwig.hmftools.isofox.fusion.cohort.FusionCohort;
import com.hartwig.hmftools.isofox.novel.cohort.AltSjCohortAnalyser;
import com.hartwig.hmftools.isofox.novel.cohort.AltSjCohortMatrix;
import com.hartwig.hmftools.isofox.novel.cohort.RecurrentVariantFinder;
import com.hartwig.hmftools.isofox.novel.cohort.SpliceSiteCache;
import com.hartwig.hmftools.isofox.novel.cohort.SpliceVariantMatcher;
import com.hartwig.hmftools.isofox.unmapped.UmrCohortAnalyser;

import org.jetbrains.annotations.NotNull;

public class CohortAnalyser
{
    private final CohortConfig mConfig;
    private final ConfigBuilder mCmdLineArgs;

    public CohortAnalyser(final ConfigBuilder configBuilder)
    {
        mCmdLineArgs = configBuilder;
        mConfig = new CohortConfig(configBuilder);
    }

    public boolean run()
    {
        if(!mConfig.SampleData.isValid())
            return false;

        for(AnalysisType type : mConfig.AnalysisTypes)
        {
            switch(type)
            {
                case SUMMARY:
                    loadSummaryData();
                    break;

                case ALT_SPLICE_JUNCTION:
                {
                    AltSjCohortAnalyser altSjCohort = new AltSjCohortAnalyser(mConfig, mCmdLineArgs);
                    altSjCohort.processAltSpliceJunctions();
                    break;
                }

                case ALT_SPLICE_JUNCTION_MATRIX:
                {
                    AltSjCohortMatrix altSjMatrix = new AltSjCohortMatrix(mConfig, mCmdLineArgs);
                    altSjMatrix.processAltSpliceJunctions();
                    break;
                }

                case SPLICE_VARIANT_MATCHING:
                {
                    SpliceVariantMatcher spliceVariantMatcher = new SpliceVariantMatcher(mConfig, mCmdLineArgs);
                    spliceVariantMatcher.processAltSpliceJunctions();
                    break;
                }

                case UNMAPPED_READS:
                {
                    UmrCohortAnalyser umrAnalyser = new UmrCohortAnalyser(mConfig, mCmdLineArgs);
                    umrAnalyser.processSampleFiles();
                    break;
                }

                case FUSION:
                {
                    FusionCohort fusionCohort = new FusionCohort(mConfig, mCmdLineArgs);
                    fusionCohort.processFusionFiles();
                    break;
                }

                case EXPRESSION_DISTRIBUTION:
                {
                    ExpressionCohortDistribution expressionDistibution = new ExpressionCohortDistribution(mConfig, mCmdLineArgs);
                    expressionDistibution.produceCohortData();
                    break;
                }

                case SPLICE_SITE_PERCENTILES:
                {
                    SpliceSiteCache spliceSiteCache = new SpliceSiteCache(mConfig, mCmdLineArgs);
                    spliceSiteCache.createPercentiles();
                    break;
                }

                case RECURRENT_SPLICE_VARIANTS:
                {
                    RecurrentVariantFinder recurrentVariantFinder = new RecurrentVariantFinder(mConfig, mCmdLineArgs);
                    recurrentVariantFinder.processCohortSomaticVariants();
                    break;
                }

                case EXTERNAL_EXPRESSION_COMPARE:
                {
                    ExternalExpressionCompare expExpressionCompare = new ExternalExpressionCompare(mConfig);
                    expExpressionCompare.processSampleTranscriptFiles();
                    break;
                }

                case GENE_EXPRESSION_COMPARE:
                {
                    ExpressionCohortCompare expCompare = new ExpressionCohortCompare(mConfig);
                    expCompare.runAnalysis();
                    break;
                }

                case GENE_EXPRESSION_MATRIX:
                case TRANSCRIPT_EXPRESSION_MATRIX:
                {
                    ExpressionMatrix expMatrix = new ExpressionMatrix(mConfig, type);
                    expMatrix.processSamples();
                    break;
                }

                case PANEL_TPM_NORMALISATION:
                {
                    GeneratePanelNormalisation generatePanelNormalisation = new GeneratePanelNormalisation(mConfig, mCmdLineArgs);
                    generatePanelNormalisation.processSamples();
                    break;
                }

                default:
                    break;
            }
        }

        ISF_LOGGER.info("Isofox cohort analyser complete");

        return true;
    }

    private void loadSummaryData()
    {
        final List<Path> filenames = Lists.newArrayList();

        if(!formSampleFilenames(mConfig, SUMMARY, filenames))
            return;

        try
        {
            // for now write out a single consolidated file
            final List<RnaStatistics> summaryStats = filenames.stream().map(x -> loadFile(x))
                    .filter(x -> x != null)
                    .collect(Collectors.toList());

            if(summaryStats.size() != mConfig.SampleData.SampleIds.size())
                return;

            final String outputFileName = mConfig.formCohortFilename(SUMMARY_FILE_ID);
            final BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write(RnaStatisticFile.header());
            writer.newLine();

            for(int i = 0; i < summaryStats.size(); ++i)
            {
                final RnaStatistics stats = summaryStats.get(i);
                writer.write(RnaStatisticFile.writeLine(mConfig.SampleData.SampleIds.get(i), stats));
                writer.newLine();
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write cohort summary data file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        CohortConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        CohortAnalyser cohortAnalyser = new CohortAnalyser(configBuilder);

        if(!cohortAnalyser.run())
        {
            ISF_LOGGER.info("Isofox cohort analyser failed");
            System.exit(1);
        }
    }
}
