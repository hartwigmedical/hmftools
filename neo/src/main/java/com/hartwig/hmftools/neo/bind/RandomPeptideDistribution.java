package com.hartwig.hmftools.neo.bind;

import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.DELIM;
import static com.hartwig.hmftools.neo.bind.BinderConfig.formFilename;
import static com.hartwig.hmftools.neo.bind.GlobalWeights.GLOBAL_COUNTS;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sigs.VectorUtils;
import com.hartwig.hmftools.common.utils.Doubles;

public class RandomPeptideDistribution
{
    private final RandomPeptideConfig mConfig;
    private boolean mDataLoaded;

    private final Map<String,Map<Integer,List<ScoreDistributionData>>> mAlleleScoresMap; // allele to peptide length to distribution
    private final List<double[]> mDiscreteScoreData;

    private static final int SCORE_SIZE = 0;
    private static final int SCORE_BRACKET = 1;

    public static final double INVALID_RANK = -1;

    public RandomPeptideDistribution(final RandomPeptideConfig config)
    {
        mConfig = config;

        mAlleleScoresMap = Maps.newHashMap();
        mDataLoaded = false;

        mDiscreteScoreData = Lists.newArrayList();
        mDiscreteScoreData.add(new double[] {0.0001, 0.01});
        mDiscreteScoreData.add(new double[] {0.001, 0.1});
        mDiscreteScoreData.add(new double[] {0.01, 0.25});
        mDiscreteScoreData.add(new double[] {0.05, 1.0});
    }

    public boolean loadData()
    {
        mDataLoaded = loadDistribution();
        return mDataLoaded;
    }

    public boolean hasData() { return mDataLoaded; }

    public Map<String,Map<Integer,List<ScoreDistributionData>>> getAlleleScoresMap() { return mAlleleScoresMap; }

    public double getScoreRank(final String allele, final int peptideLength, double score)
    {
        Map<Integer,List<ScoreDistributionData>> peptideLengthMap = mAlleleScoresMap.get(allele);

        if(peptideLengthMap == null)
            return INVALID_RANK;

        List<ScoreDistributionData> scores = peptideLengthMap.get(peptideLength);

        if(scores == null || scores.size() < 2)
            return INVALID_RANK;

        boolean isAscending = scores.get(0).Score < scores.get(1).Score;

        if((isAscending && score < scores.get(0).Score) || (!isAscending && score > scores.get(0).Score))
            return 0; // zero-th percentile if the score is better than any in the random distribution

        for(int i = 0; i < scores.size(); ++i)
        {
            ScoreDistributionData scoreData = scores.get(i);

            if(Doubles.equal(score, scoreData.Score))
                return scoreData.ScoreBucket;

            ScoreDistributionData nextScoreData = i < scores.size() - 1 ? scores.get(i + 1) : null;

            if(Doubles.equal(score, nextScoreData.Score))
                return nextScoreData.ScoreBucket;

            if((isAscending && score > scoreData.Score) || (!isAscending && score < scoreData.Score))
            {
                if(nextScoreData == null)
                    break;

                if((isAscending && score < nextScoreData.Score) || (!isAscending && score > nextScoreData.Score))
                    return (scoreData.ScoreBucket + nextScoreData.ScoreBucket) * 0.5;
            }
        }

        return 1;
    }

    private List<String> loadRandomPeptides()
    {
        List<String> refRandomPeptides = Lists.newArrayList();

        if(mConfig.RandomPeptidesFile == null)
        {
            NE_LOGGER.error("missing random peptides file");
            return refRandomPeptides;
        }

        try
        {
            refRandomPeptides.addAll(Files.readAllLines(new File(mConfig.RandomPeptidesFile).toPath()));
            refRandomPeptides.remove(0);

            NE_LOGGER.info("loaded {} random peptides", refRandomPeptides.size());
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to load random peptide file: {}", e.toString());
        }

        return refRandomPeptides;
    }

    public void buildDistribution(final Map<String,Map<Integer,BindScoreMatrix>> alleleBindMatrixMap)
    {
        List<String> refRandomPeptides = loadRandomPeptides();

        if(refRandomPeptides.isEmpty())
            return;

        // score each against each allele and build up a percentiles for each
        int alleleCount = 0;

        for(Map.Entry<String,Map<Integer,BindScoreMatrix>> alleleEntry : alleleBindMatrixMap.entrySet())
        {
            final String allele = alleleEntry.getKey();

            if(!allele.equals(GLOBAL_COUNTS) && !mConfig.AllelesToWrite.isEmpty() && !mConfig.AllelesToWrite.contains(allele))
                continue;

            NE_LOGGER.debug("building distribution for allele({})", allele);

            final Map<Integer,BindScoreMatrix> peptideLengthMatrixMap = alleleEntry.getValue();

            Map<Integer,List<ScoreDistributionData>> peptideLengthMap = Maps.newHashMap();
            mAlleleScoresMap.put(allele, peptideLengthMap);

            for(BindScoreMatrix matrix : peptideLengthMatrixMap.values())
            {

                List<Double> peptideScores = Lists.newArrayListWithExpectedSize(refRandomPeptides.size());

                int count = 0;

                for(String peptide : refRandomPeptides)
                {
                    if(peptide.length() > matrix.PeptideLength)
                        peptide = peptide.substring(0, matrix.PeptideLength);

                    double score = matrix.calcScore(peptide);
                    VectorUtils.optimisedAdd(peptideScores, score, false);

                    ++count;

                    if(count > 0 && (count % 500000) == 0)
                    {
                        NE_LOGGER.debug("added {} sorted random peptide scores", count);
                    }
                }

                List<ScoreDistributionData> scoresDistributions = generateDistribution(matrix.Allele, matrix.PeptideLength, peptideScores);
                peptideLengthMap.put(matrix.PeptideLength, scoresDistributions);
            }

            ++alleleCount;

            if(alleleCount > 0 && (alleleCount % 10) == 0)
            {
                NE_LOGGER.info("generated distributions for {} alleles", alleleCount);
            }
        }

        if(mConfig.WriteRandomDistribution)
            writeDistribution();
    }

    public List<ScoreDistributionData> generateDistribution(final String allele, final int peptideLength, final List<Double> peptideScores)
    {
        // write the distribution as 0.0001 up to 0.01, 0.001 up to 0.01, then 0.01 up to 100%
        int totalScores = peptideScores.size();
        int discreteIndex = 0;
        double currentSize = mDiscreteScoreData.get(discreteIndex)[SCORE_SIZE];
        double currentBracket = mDiscreteScoreData.get(discreteIndex)[SCORE_BRACKET];
        int requiredScores = (int) round(totalScores * currentSize);

        double scoreTotal = 0;
        int currentScoreCount = 0;
        double currentSizeTotal = 0;
        int cumulativeScores = 0;

        List<ScoreDistributionData> scoresDistributions = Lists.newArrayList();

        for(Double score : peptideScores)
        {
            scoreTotal += score;
            ++currentScoreCount;
            ++cumulativeScores;

            if(currentScoreCount >= requiredScores)
            {
                double avgScore = scoreTotal / currentScoreCount;

                scoresDistributions.add(new ScoreDistributionData(
                        allele, peptideLength, currentSizeTotal, avgScore, currentScoreCount, cumulativeScores));

                scoreTotal = 0;
                currentScoreCount = 0;
                currentSizeTotal += currentSize;

                if(Doubles.equal(currentSizeTotal, currentBracket))
                {
                    ++discreteIndex;

                    if(discreteIndex < mDiscreteScoreData.size())
                    {
                        currentSize = mDiscreteScoreData.get(discreteIndex)[SCORE_SIZE];
                        currentBracket = mDiscreteScoreData.get(discreteIndex)[SCORE_BRACKET];
                        requiredScores = (int) round(totalScores * currentSize);
                    }
                }
            }
        }

        return scoresDistributions;
    }

    public static BufferedWriter initialiseWriter(final String filename)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Allele,PeptideLength,ScoreBucket,Score,BucketCount,CumulativeCount");
            // writer.write(",CurrentScores,CurrentSize,CurrentBracket,CurrentSizeTotal");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write random peptide file: {}", e.toString());
            return null;
        }
    }

    public static void writeDistribution(final BufferedWriter writer, final List<ScoreDistributionData> scoreDistributionData)
    {
        try
        {
            for(ScoreDistributionData scoreDist : scoreDistributionData)
            {
                writer.write(String.format("%s,%d,%f,%.4f,%d,%d",
                        scoreDist.Allele, scoreDist.PeptideLength, scoreDist.ScoreBucket, scoreDist.Score,
                        scoreDist.BucketCount, scoreDist.CumulativeCount));

                //writer.write(String.format(",%d,%f,%f,%f", currentScores, currentSize, currentBracket, currentSizeTotal));
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write random peptide file: {}", e.toString());
            return;
        }
    }

    private void writeDistribution()
    {
        NE_LOGGER.info("writing random peptide scoring distribution for {} alleles",
                !mConfig.AllelesToWrite.isEmpty() ? mConfig.AllelesToWrite.size() : "all");

        final String distFilename = formFilename("random_peptide_dist", mConfig.OutputDir, mConfig.OutputId);

        try
        {
            BufferedWriter writer = initialiseWriter(distFilename);

            // score each against each allele and build up a percentiles for each
            for(Map.Entry<String,Map<Integer,List<ScoreDistributionData>>> alleleEntry : mAlleleScoresMap.entrySet())
            {
                String allele = alleleEntry.getKey();

                final Map<Integer,List<ScoreDistributionData>> pepLenDistMap = alleleEntry.getValue();

                for(List<ScoreDistributionData> pepLenScoreDist : pepLenDistMap.values())
                {
                    writeDistribution(writer, pepLenScoreDist);
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write random peptide file: {}", e.toString());
        }
    }

    private boolean loadDistribution()
    {
        if(mConfig.RandomPeptideDistributionFile == null)
            return false;

        try
        {
            final List<String> lines = Files.readAllLines(new File(mConfig.RandomPeptideDistributionFile).toPath());

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIM);
            lines.remove(0);

            for(String line : lines)
            {
                ScoreDistributionData data = ScoreDistributionData.fromCsv(line, fieldsIndexMap);

                Map<Integer,List<ScoreDistributionData>> peptideLengthMap = mAlleleScoresMap.get(data.Allele);

                if(peptideLengthMap == null)
                {
                    peptideLengthMap = Maps.newHashMap();
                    mAlleleScoresMap.put(data.Allele, peptideLengthMap);
                }

                List<ScoreDistributionData> dataList = peptideLengthMap.get(data.PeptideLength);

                if(dataList == null)
                {
                    dataList = Lists.newArrayList();
                    peptideLengthMap.put(data.PeptideLength, dataList);
                }

                dataList.add(data);
            }

            NE_LOGGER.info("loaded {} alleles with random peptide score distribution data", mAlleleScoresMap.size());
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to read random distribution file: {}", e.toString());
            return false;
        }

        return true;
    }
}
