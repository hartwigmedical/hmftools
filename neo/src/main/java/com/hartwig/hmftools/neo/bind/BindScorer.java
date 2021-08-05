package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindData.loadBindData;
import static com.hartwig.hmftools.neo.bind.PeptideWriteType.LIKELY_INCORRECT;
import static com.hartwig.hmftools.neo.bind.PeptideWriteType.TRAINING;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.stats.AucCalc;
import com.hartwig.hmftools.common.stats.AucData;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.jetbrains.annotations.NotNull;

public class BindScorer
{
    private final BinderConfig mConfig;

    private final Map<String,Map<Integer,List<BindData>>> mAllelePeptideData;
    private final Map<String,Map<Integer,BindScoreMatrix>> mAlleleBindMatrices;

    private final RandomPeptideDistribution mRandomDistribution;

    private final List<Double> mScoreRankPercentileBuckets;
    private final Map<String,Map<Integer,int[]>> mAlleleScoreRankPercentileCounts;


    public BindScorer(final BinderConfig config)
    {
        mConfig = config;

        mAllelePeptideData = Maps.newHashMap();
        mAlleleBindMatrices = Maps.newHashMap();
        mRandomDistribution = new RandomPeptideDistribution(config.RandomPeptides);

        // only used for training
        mScoreRankPercentileBuckets = null;
        mAlleleScoreRankPercentileCounts = null;
    }

    public BindScorer(
            final BinderConfig config, final Map<String,Map<Integer,List<BindData>>> allelePeptideData,
            final Map<String,Map<Integer,BindScoreMatrix>> alleleBindMatrices, final RandomPeptideDistribution randomDistribution)
    {
        mConfig = config;

        mAllelePeptideData = allelePeptideData;
        mAlleleBindMatrices = alleleBindMatrices;
        mRandomDistribution = randomDistribution;

        mAlleleScoreRankPercentileCounts = Maps.newHashMap();
        mScoreRankPercentileBuckets = Lists.newArrayListWithExpectedSize(16);

        double rankPercBucket = 0.00005;
        while(rankPercBucket < 1)
        {
            mScoreRankPercentileBuckets.add(rankPercBucket);
            rankPercBucket *= 2;
        }
    }

    public void run()
    {
        NE_LOGGER.info("running BindScorer");

        if(!loadData())
        {
            NE_LOGGER.error("invalid data");
            return;
        }

        runScoring();

        NE_LOGGER.info("Bind scoring complete");
    }

    public void runScoring()
    {
        NE_LOGGER.info("writing allele summaries");

        BufferedWriter alleleWriter = initAlleleSummaryWriter();

        // rank both the training and random peptide data using the newly created data and the random peptide distributions
        for(Map.Entry<String,Map<Integer,List<BindData>>> alleleEntry : mAllelePeptideData.entrySet())
        {
            final String allele = alleleEntry.getKey();
            final Map<Integer,List<BindData>> pepLenBindDataMap = alleleEntry.getValue();

            Map<Integer,BindScoreMatrix> pepLenMatrixMap = mAlleleBindMatrices.get(allele);

            if(pepLenMatrixMap == null)
            {
                NE_LOGGER.warn("allele({}) has no matrix scoring data", allele);
                continue;
            }

            for(Map.Entry<Integer,List<BindData>> pepLenEntry : pepLenBindDataMap.entrySet())
            {
                final List<BindData> bindDataList = pepLenEntry.getValue();
                for(BindData bindData : bindDataList)
                {
                    BindScoreMatrix matrix = pepLenMatrixMap.get(bindData.peptideLength());

                    if(matrix == null)
                        continue;

                    double score = matrix.calcScore(bindData.Peptide);
                    bindData.setScoreData(score, mRandomDistribution.getScoreRank(allele, bindData.peptideLength(), score));
                }
            }

            writeAlleleSummary(alleleWriter, allele);
        }

        closeBufferedWriter(alleleWriter);

        NE_LOGGER.info("writing peptide scores");
        writePeptideScores();
    }

    public void buildScoreDistribution()
    {
        if(mScoreRankPercentileBuckets == null)
            return;

        for(Map.Entry<String,Map<Integer,List<BindData>>> alleleEntry : mAllelePeptideData.entrySet())
        {
            final String allele = alleleEntry.getKey();
            final Map<Integer,List<BindData>> pepLenBindDataMap = alleleEntry.getValue();

            Map<Integer,int[]> pepLenScoreDistMap = Maps.newHashMap();
            mAlleleScoreRankPercentileCounts.put(allele, pepLenScoreDistMap);

            for(Map.Entry<Integer,List<BindData>> pepLenEntry : pepLenBindDataMap.entrySet())
            {
                int peptideLength = pepLenEntry.getKey();

                int[] distCounts = new int[mScoreRankPercentileBuckets.size()];
                pepLenScoreDistMap.put(peptideLength, distCounts);

                final List<BindData> bindDataList = pepLenEntry.getValue();
                for(BindData bindData : bindDataList)
                {
                    for(int scoreBucket = 0; scoreBucket < mScoreRankPercentileBuckets.size(); ++scoreBucket)
                    {
                        if(bindData.score() < mScoreRankPercentileBuckets.get(0) || scoreBucket == mScoreRankPercentileBuckets.size() - 1)
                        {
                            ++distCounts[scoreBucket];
                            break;
                        }
                    }
                }
            }
        }
    }

    private BufferedWriter initAlleleSummaryWriter()
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(mConfig.formFilename("allele_summary"), false);
            writer.write("Allele,PeptideLength,Source,TrainingCount,RandomCount,AUC");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to init allele summary writer: {}", e.toString());
            return null;
        }
    }

    private void writeAlleleSummary(final BufferedWriter writer, final String allele)
    {
        if(writer == null)
            return;

        try
        {
            final Map<Integer,List<BindData>> pepLenBindDataMap = mAllelePeptideData.get(allele);

            int errors = 0;

            // internal model assesses AUC per peptide length
            // external models across an entire allele

            List<AucData> alleleAucMcfData = Lists.newArrayList();
            int trainingCountMcf = 0;
            int randomCountMcf = 0;

            for(Map.Entry<Integer,List<BindData>> pepLenEntry : pepLenBindDataMap.entrySet())
            {
                List<AucData> alleleAucData = Lists.newArrayList();
                List<AucData> alleleAucRankData = Lists.newArrayList();

                int trainingCount = 0;
                int randomCount = 0;

                for(BindData bindData : pepLenEntry.getValue())
                {
                    if(bindData.isTraining())
                        ++trainingCount;
                    else
                        ++randomCount;

                    boolean isPositive = bindData.Affinity < mConfig.Constants.BindingAffinityHigh;

                    alleleAucData.add(new AucData(isPositive, bindData.score(), false));
                    alleleAucRankData.add(new AucData(isPositive, bindData.rankPercentile(), true));

                    if(bindData.hasPredictionData())
                        alleleAucMcfData.add(new AucData(isPositive, bindData.affinityPercentile(), true));
                }

                double aucPerc = AucCalc.calcPercentilesAuc(alleleAucRankData, Level.TRACE);

                NE_LOGGER.debug(String.format("allele(%s) peptideLength(%d) peptides(train=%d, rand=%d) AUC(%.4f)",
                        allele, pepLenEntry.getKey(), trainingCount, randomCount, aucPerc));

                writer.write(String.format("%s,%d,%s,%d,%d,%.4f",
                        allele, pepLenEntry.getKey(), "MODEL", trainingCount, randomCount, aucPerc));
                writer.newLine();

                trainingCountMcf += trainingCount;
                randomCountMcf += randomCount;
            }

            if(!alleleAucMcfData.isEmpty())
            {
                double aucMcf = AucCalc.calcPercentilesAuc(alleleAucMcfData, Level.TRACE);

                NE_LOGGER.debug(String.format("allele(%s) McFlurry peptides(train=%d, rand=%d) AUC(%.4f)",
                        allele, trainingCountMcf, randomCountMcf, aucMcf));

                writer.write(String.format("%s,%d,%s,%d,%d,%.4f",
                        allele, 0, "MCF", trainingCountMcf, randomCountMcf, aucMcf));
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write allele summary data: {}", e.toString());
        }
    }

    private void writePeptideScores()
    {
        if(mConfig.WritePeptideType == PeptideWriteType.NONE)
            return;

        try
        {
            BufferedWriter writer = createBufferedWriter(mConfig.formFilename("peptide_scores"), false);
            writer.write("Allele,Peptide,Source,Score,Rank,Affinity,PredictedAffinity,AffinityPerc,PresentationPerc");
            writer.newLine();

            for(Map.Entry<String,Map<Integer,List<BindData>>> alleleEntry : mAllelePeptideData.entrySet())
            {
                final String allele = alleleEntry.getKey();
                final Map<Integer,List<BindData>> pepLenBindDataMap = alleleEntry.getValue();

                for(List<BindData> bindDataList : pepLenBindDataMap.values())
                {
                    for(BindData bindData : bindDataList)
                    {
                        if(mConfig.WritePeptideType == LIKELY_INCORRECT)
                        {
                            boolean expectBinds = bindData.Affinity < mConfig.Constants.BindingAffinityHigh;
                            boolean inTopPerc = bindData.rankPercentile() < 0.02;
                            if(expectBinds == inTopPerc)
                                continue;
                        }
                        else if(mConfig.WritePeptideType == TRAINING && !bindData.isTraining())
                        {
                            continue;
                        }

                        writer.write(String.format("%s,%s,%s,%.4f,%.4f,%.2f,%.2f,%.4f,%.4f",
                                allele, bindData.Peptide, bindData.Source,
                                bindData.score(), bindData.rankPercentile(), bindData.Affinity,
                                bindData.predictedAffinity(), bindData.affinityPercentile(), bindData.presentationPercentile()));
                        writer.newLine();
                    }
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write peptide scores file: {}", e.toString());
        }
    }

    private boolean loadData()
    {
        if(!loadBindData(
                mConfig.TrainingDataFile, true, mConfig.SpecificAlleles, mConfig.SpecificPeptideLengths, mAllelePeptideData))
        {
            return false;
        }

        List<BindScoreMatrix> matrixList = BindScoreMatrix.loadFromCsv(mConfig.BindMatrixFile);

        for(BindScoreMatrix matrix : matrixList)
        {
            Map<Integer,BindScoreMatrix> pepLenMap = mAlleleBindMatrices.get(matrix.Allele);

            if(pepLenMap == null)
            {
                pepLenMap = Maps.newHashMap();
                mAlleleBindMatrices.put(matrix.Allele, pepLenMap);
            }

            pepLenMap.put(matrix.PeptideLength, matrix);
        }

        if(!mRandomDistribution.loadData())
            return false;

        return true;
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        BinderConfig.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        BinderConfig config = new BinderConfig(cmd);

        BindScorer bindScorer = new BindScorer(config);
        bindScorer.run();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
