package com.hartwig.hmftools.isofox.cohort;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConfig.LOG_DEBUG;
import static com.hartwig.hmftools.isofox.cohort.CohortAnalysisType.SUMMARY;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.isValid;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.SUMMARY_FILE;
import static com.hartwig.hmftools.isofox.results.SummaryStats.loadFile;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.isofox.expression.cohort.ExpressionCohortCompare;
import com.hartwig.hmftools.isofox.expression.cohort.ExpressionMatrix;
import com.hartwig.hmftools.isofox.expression.cohort.ExternalExpressionCompare;
import com.hartwig.hmftools.isofox.expression.cohort.GeneExpressionDistribution;
import com.hartwig.hmftools.isofox.expression.cohort.SampleGenePercentiles;
import com.hartwig.hmftools.isofox.expression.cohort.TransExpressionDistribution;
import com.hartwig.hmftools.isofox.fusion.cohort.FusionCohort;
import com.hartwig.hmftools.isofox.novel.cohort.AltSjCohortAnalyser;
import com.hartwig.hmftools.isofox.novel.cohort.SpliceVariantMatcher;
import com.hartwig.hmftools.isofox.results.SummaryStats;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class CohortAnalyser
{
    private final CohortConfig mConfig;
    private final CommandLine mCmdLineArgs;

    public CohortAnalyser(final CommandLine cmdLineArgs)
    {
        mCmdLineArgs = cmdLineArgs;
        mConfig = new CohortConfig(cmdLineArgs);
    }

    public boolean load()
    {
        if(!mConfig.SampleData.isValid())
            return false;

        for(CohortAnalysisType type : mConfig.AnalysisTypes)
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

                case SPLICE_VARIANT_MATCHING:
                {
                    SpliceVariantMatcher spliceVariantMatcher = new SpliceVariantMatcher(mConfig, mCmdLineArgs);
                    spliceVariantMatcher.processAltSpliceJunctions();
                    break;
                }

                case FUSION:
                {
                    FusionCohort fusionCohort = new FusionCohort(mConfig, mCmdLineArgs);
                    fusionCohort.processFusionFiles();
                    break;
                }

                case TRANSCRIPT_DISTRIBUTION:
                {
                    TransExpressionDistribution transExpDist = new TransExpressionDistribution(mConfig);
                    transExpDist.processSampleTranscriptFiles();
                    break;
                }

                case GENE_DISTRIBUTION:
                {
                    GeneExpressionDistribution geneExpDistribution = new GeneExpressionDistribution(mConfig);
                    geneExpDistribution.processGenes();
                    break;
                }

                case SAMPLE_GENE_PERCENTILES:
                {
                    SampleGenePercentiles sampleGenePerc = new SampleGenePercentiles(mConfig);
                    sampleGenePerc.processSampleFiles();
                    break;
                }

                case SAMPLE_ROUTINES:
                {
                    SampleRoutines sampleRoutines = new SampleRoutines(mConfig);
                    sampleRoutines.convertMutationCohortFile(mConfig.SampleMutationsFile);
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
                    expCompare.processSamples();
                    break;
                }

                case GENE_EXPRESSION_MATRIX:
                case TRANSCRIPT_EXPRESSION_MATRIX:
                {
                    ExpressionMatrix expMatrix = new ExpressionMatrix(mConfig, type);
                    expMatrix.processSamples();
                    break;
                }

                default:
                    break;
            }
        }

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
            final List<SummaryStats> summaryStats = filenames.stream().map(x -> loadFile(x))
                    .filter(x -> x != null)
                    .collect(Collectors.toList());

            if(summaryStats.size() != mConfig.SampleData.SampleIds.size())
                return;

            final String outputFileName = mConfig.formCohortFilename(SUMMARY_FILE);
            final BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write(SummaryStats.csvHeader());
            writer.newLine();

            for(int i = 0; i < summaryStats.size(); ++i)
            {
                final SummaryStats stats = summaryStats.get(i);
                writer.write(stats.toCsv(mConfig.SampleData.SampleIds.get(i)));
                writer.newLine();
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write cohort summary data file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = CohortConfig.createCmdLineOptions();
        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        if(!isValid(cmd))
        {
            ISF_LOGGER.error("missing or invalid config options");
            return;
        }

        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        CohortAnalyser cohortAnalyser = new CohortAnalyser(cmd);

        if(!cohortAnalyser.load())
        {
            ISF_LOGGER.info("Isofox cohort analyser failed");
            return;
        }

        ISF_LOGGER.info("Isofox cohort analyser complete");
    }

}
