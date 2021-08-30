package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindData.loadBindData;

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
    private final FlankScores mFlankScores;

    private final RandomPeptideDistribution mRandomDistribution;
    private final BindingLikelihood mBindingLikelihood;

    public BindScorer(final BinderConfig config)
    {
        mConfig = config;

        mAllelePeptideData = Maps.newHashMap();
        mAlleleBindMatrices = Maps.newHashMap();
        mRandomDistribution = new RandomPeptideDistribution(config.RandomPeptides);

        mBindingLikelihood = new BindingLikelihood();
        mFlankScores = new FlankScores();
    }

    public BindScorer(
            final BinderConfig config, final Map<String,Map<Integer,List<BindData>>> allelePeptideData,
            final Map<String,Map<Integer,BindScoreMatrix>> alleleBindMatrices, final RandomPeptideDistribution randomDistribution,
            final FlankScores flankScores)
    {
        mConfig = config;

        mAllelePeptideData = allelePeptideData;
        mAlleleBindMatrices = alleleBindMatrices;
        mRandomDistribution = randomDistribution;
        mBindingLikelihood = null;
        mFlankScores = flankScores;
    }

    public void run()
    {
        NE_LOGGER.info("running BindScorer");

        if(!loadData())
        {
            NE_LOGGER.error("invalid reference data");
            return;
        }

        if(!loadBindData(
                mConfig.ValidationDataFile, true, Lists.newArrayList(), Lists.newArrayList(), mAllelePeptideData))
        {
            NE_LOGGER.error("invalid validation data");
            return;
        }

        runScoring();

        NE_LOGGER.info("scoring complete");
    }

    public static final double INVALID_CALC = -1;

    public double calcScore(final String allele, final String peptide)
    {
        Map<Integer,BindScoreMatrix> pepLenMatrixMap = mAlleleBindMatrices.get(allele);
        if(pepLenMatrixMap == null)
            return INVALID_CALC;

        BindScoreMatrix matrix = pepLenMatrixMap.get(peptide.length());

        if(matrix == null)
            return INVALID_CALC;

        return matrix.calcScore(peptide);
    }

    public double calcScoreRank(final String allele, final String peptide, double score)
    {
        return mRandomDistribution.getScoreRank(allele, peptide.length(), score);
    }

    public double calcLikelihood(final String allele, final String peptide, double ranks)
    {
        if(mBindingLikelihood == null || !mBindingLikelihood.hasData())
            return INVALID_CALC;

        return mBindingLikelihood.getBindingLikelihood(allele, peptide, ranks);
    }

    public void runScoring()
    {
        NE_LOGGER.info("running scoring");

        BufferedWriter alleleWriter = mConfig.WriteSummaryData ? initAlleleSummaryWriter() : null;

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

            TprCalc alleleTprCalc = new TprCalc();

            for(Map.Entry<Integer,List<BindData>> pepLenEntry : pepLenBindDataMap.entrySet())
            {
                final List<BindData> bindDataList = pepLenEntry.getValue();

                for(BindData bindData : bindDataList)
                {
                    BindScoreMatrix matrix = pepLenMatrixMap.get(bindData.peptideLength());

                    if(matrix == null)
                        continue;

                    calcScoreData(bindData, matrix, mFlankScores, mRandomDistribution, mBindingLikelihood);

                    alleleTprCalc.addRank(bindData.likelihoodRank());
                }
            }

            if(mConfig.WriteSummaryData)
                writeAlleleSummary(alleleWriter, allele, alleleTprCalc);
        }

        closeBufferedWriter(alleleWriter);

        writePeptideScores();
    }

    public static double calcScore(
            final BindScoreMatrix matrix, final FlankScores flankScores, final String peptide, final String upFlank, final String downFlank)
    {
        double score = matrix.calcScore(peptide);

        if(flankScores.hasData())
        {
            double flankScore = flankScores.calcScore(upFlank, downFlank);
            score += flankScore;
        }

        return score;
    }

    public static void calcScoreData(
            final BindData bindData, final BindScoreMatrix matrix, final FlankScores flankScores,
            final RandomPeptideDistribution randomDistribution, final BindingLikelihood bindingLikelihood)
    {
        double score = matrix.calcScore(bindData.Peptide);

        double flankScore = 0;
        if(flankScores.hasData() && bindData.hasFlanks())
        {
            flankScore = flankScores.calcScore(bindData.UpFlank, bindData.DownFlank);
            score += flankScore;
        }

        double rankPercentile = randomDistribution.getScoreRank(bindData.Allele, bindData.peptideLength(), score);

        double likelihood = bindingLikelihood != null && bindingLikelihood.hasData() ?
                bindingLikelihood.getBindingLikelihood(bindData.Allele, bindData.Peptide, rankPercentile) : 0;

        double likelihoodRank = randomDistribution.getLikelihoodRank(bindData.Allele, likelihood);

        bindData.setScoreData(score, flankScore, rankPercentile, likelihood, likelihoodRank);
    }

    private void writePeptideScores()
    {
        if(mConfig.WritePeptideType == PeptideWriteType.NONE)
            return;

        NE_LOGGER.info("writing peptide scores");

        try
        {
            BufferedWriter writer = createBufferedWriter(mConfig.formOutputFilename("peptide_scores"), false);
            writer.write("Allele,Peptide,Source,Score,Rank,Likelihood,LikelihoodRank");

            if(mConfig.ApplyFlanks)
                writer.write(",FlankScore,UpFlank,DownFlank");

            boolean hasPredictionData = false;
            boolean hasMeasuredAffinity = false;

            if(mAllelePeptideData.values().stream().anyMatch(x -> x.values().stream().filter(y -> !y.isEmpty()).anyMatch(y -> y.get(0).hasMeasuredAffinity())))
            {
                hasMeasuredAffinity = true;
                writer.write(",MeasuredAffinity");
            }

            if(mAllelePeptideData.values().stream().anyMatch(x -> x.values().stream().filter(y -> !y.isEmpty()).anyMatch(y -> y.get(0).hasPredictionData())))
            {
                hasPredictionData = true;
                writer.write(",PredictedAffinity,AffinityPerc,PresScore,PresPerc");
            }

            writer.newLine();

            for(Map.Entry<String,Map<Integer,List<BindData>>> alleleEntry : mAllelePeptideData.entrySet())
            {
                final String allele = alleleEntry.getKey();
                final Map<Integer,List<BindData>> pepLenBindDataMap = alleleEntry.getValue();

                for(List<BindData> bindDataList : pepLenBindDataMap.values())
                {
                    for(BindData bindData : bindDataList)
                    {
                        writer.write(String.format("%s,%s,%s,%.4f,%.6f,%.6f,%.6f",
                                allele, bindData.Peptide, bindData.Source, bindData.score(),
                                bindData.rankPercentile(), bindData.likelihood(), bindData.likelihoodRank()));

                        if(mConfig.ApplyFlanks)
                        {
                            writer.write(String.format(",%.4f,%s,%s", bindData.flankScore(), bindData.UpFlank, bindData.DownFlank));
                        }

                        if(hasMeasuredAffinity)
                        {
                            writer.write(String.format(",%.2f", bindData.measuredAffinity()));
                        }

                        if(hasPredictionData)
                        {
                            writer.write(String.format(",%.2f,%.6f,%.4f,%.6f",
                                    bindData.predictedAffinity(), bindData.affinityPercentile(),
                                    bindData.presentationScore(), bindData.presentationPercentile()));
                        }

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

    private BufferedWriter initAlleleSummaryWriter()
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(mConfig.formOutputFilename("allele_summary"), false);
            writer.write("Allele,BindCount,TPR");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to init allele summary writer: {}", e.toString());
            return null;
        }
    }

    private void writeAlleleSummary(final BufferedWriter writer, final String allele, final TprCalc alleleTprCalc)
    {
        if(writer == null)
            return;

        try
        {
            writer.write(String.format("%s,%d,%.4f",
                    allele, alleleTprCalc.entryCount(), alleleTprCalc.calc()));
            writer.newLine();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write allele summary data: {}", e.toString());
        }
    }

    private void writeAlleleSummary(final BufferedWriter writer, final String allele)
    {
        if(writer == null)
            return;

        try
        {
            final Map<Integer,List<BindData>> pepLenBindDataMap = mAllelePeptideData.get(allele);

            // internal model assesses AUC per peptide length
            // external models across an entire allele

            List<AucData> alleleAucMcfData = Lists.newArrayList();
            int trainingCountMcf = 0;

            for(Map.Entry<Integer,List<BindData>> pepLenEntry : pepLenBindDataMap.entrySet())
            {
                List<AucData> alleleAucData = Lists.newArrayList();
                List<AucData> alleleAucRankData = Lists.newArrayList();

                int trainingCount = 0;

                for(BindData bindData : pepLenEntry.getValue())
                {
                    if(bindData.isTraining())
                        ++trainingCount;

                    alleleAucData.add(new AucData(true, bindData.score(), false));
                    alleleAucRankData.add(new AucData(true, bindData.rankPercentile(), true));

                    if(bindData.hasPredictionData())
                        alleleAucMcfData.add(new AucData(true, bindData.affinityPercentile(), true));
                }

                double aucPerc = AucCalc.calcPercentilesAuc(alleleAucRankData, Level.TRACE);

                NE_LOGGER.debug(String.format("allele(%s) peptideLength(%d) peptides(%d) AUC(%.4f)",
                        allele, pepLenEntry.getKey(), trainingCount, aucPerc));

                writer.write(String.format("%s,%d,%s,%d,%.4f",
                        allele, pepLenEntry.getKey(), "MODEL", trainingCount, aucPerc));
                writer.newLine();

                trainingCountMcf += trainingCount;
            }

            if(!alleleAucMcfData.isEmpty())
            {
                double aucMcf = AucCalc.calcPercentilesAuc(alleleAucMcfData, Level.TRACE);

                NE_LOGGER.debug(String.format("allele(%s) McFlurry peptides(%d) AUC(%.4f)",
                        allele, trainingCountMcf, aucMcf));
            }
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write allele summary data: {}", e.toString());
        }
    }

    public boolean loadData()
    {
        NE_LOGGER.info("loading matrix data from {}", mConfig.PosWeightsFile);

        List<BindScoreMatrix> matrixList = BindScoreMatrix.loadFromCsv(mConfig.PosWeightsFile);

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

        NE_LOGGER.info("loading random distribution data from {}", mConfig.RandomPeptides.RandomPeptideDistributionFile);

        if(!mRandomDistribution.loadData())
            return false;

        if(mConfig.BindLikelihoodFile != null && !mBindingLikelihood.loadLikelihoods(mConfig.BindLikelihoodFile))
            return false;

        if(mConfig.FlankPosWeightsFile != null && mConfig.ApplyFlanks && !mFlankScores.loadPosWeights(mConfig.FlankPosWeightsFile))
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
