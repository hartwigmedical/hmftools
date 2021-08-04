package com.hartwig.hmftools.neo.bind;

import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindData.DELIM;

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
    private final BinderConfig mConfig;
    private boolean mDataLoaded;

    private final Map<String,Map<Integer,List<ScoreDistributionData>>> mAlleleScoresMap; // allele to peptide length to distribution

    private static final List<Double> SCORE_DISCRETE_SIZES = Lists.newArrayList(0.0001, 0.001, 0.01);
    private static final List<Double> SCORE_DISCRETE_BRACKETS = Lists.newArrayList(0.01, 0.1, 1.0);

    public static final double INVALID_RANK = -1;

    public RandomPeptideDistribution(final BinderConfig config)
    {
        mConfig = config;

        mAlleleScoresMap = Maps.newHashMap();
        mDataLoaded = false;
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

        if(scores == null || scores.isEmpty())
            return INVALID_RANK;

        if(score > scores.get(0).Score)
            return 0; // zero-th percentile if the score is better than any in the random distribution

        for(int i = 0; i < scores.size(); ++i)
        {
            ScoreDistributionData scoreData = scores.get(i);
            ScoreDistributionData nextScoreData = i < scores.size() - 1 ? scores.get(i + 1) : null;

            if(score < scoreData.Score)
            {
                if(nextScoreData == null)
                    break;

                if(score > nextScoreData.Score)
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
        for(Map.Entry<String,Map<Integer,BindScoreMatrix>> alleleEntry : alleleBindMatrixMap.entrySet())
        {
            final String allele = alleleEntry.getKey();
            final Map<Integer,BindScoreMatrix> peptideLengthMatrixMap = alleleEntry.getValue();

            Map<Integer,List<ScoreDistributionData>> peptideLengthMap = Maps.newHashMap();
            mAlleleScoresMap.put(allele, peptideLengthMap);

            for(BindScoreMatrix matrix : peptideLengthMatrixMap.values())
            {
                NE_LOGGER.info("building distribution for allele({}) peptideLength({})", matrix.Allele, matrix.PeptideLength);

                int peptideCount = refRandomPeptides.size();

                List<Double> peptideScores = Lists.newArrayList();

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

                // write the distribution as 0.0001 up to 0.01, 0.001 up to 0.01, then 0.01 up to 100%
                int discreteIndex = 0;
                double currentSize = SCORE_DISCRETE_SIZES.get(discreteIndex);
                double currentBracket = SCORE_DISCRETE_BRACKETS.get(discreteIndex);
                int requiredScores = (int) round(peptideCount * currentSize);

                double scoreTotal = 0;
                int currentScoreCount = 0;
                double currentSizeTotal = 0;
                int cumulativeScores = 0;

                List<ScoreDistributionData> scoresDistributions = Lists.newArrayList();
                peptideLengthMap.put(matrix.PeptideLength, scoresDistributions);

                for(Double score : peptideScores)
                {
                    scoreTotal += score;
                    ++currentScoreCount;
                    ++cumulativeScores;

                    if(currentScoreCount >= requiredScores)
                    {
                        double avgScore = scoreTotal / currentScoreCount;

                        scoresDistributions.add(new ScoreDistributionData(
                                matrix.Allele, matrix.PeptideLength, currentSizeTotal, avgScore, currentScoreCount, cumulativeScores));

                        scoreTotal = 0;
                        currentScoreCount = 0;
                        currentSizeTotal += currentSize;

                        if(Doubles.equal(currentSizeTotal, currentBracket))
                        {
                            ++discreteIndex;

                            if(discreteIndex < SCORE_DISCRETE_SIZES.size())
                            {
                                currentSize = SCORE_DISCRETE_SIZES.get(discreteIndex);
                                currentBracket = SCORE_DISCRETE_BRACKETS.get(discreteIndex);
                                requiredScores = (int) round(peptideCount * currentSize);
                            }
                        }
                    }
                }
            }
        }

        if(mConfig.WriteRandomDistribution)
            writeDistribution();
    }

    private void writeDistribution()
    {
        NE_LOGGER.info("writing random peptide scoring distribution");

        final String distFilename = mConfig.formFilename("random_peptide_dist");

        try
        {
            BufferedWriter writer = createBufferedWriter(distFilename, false);

            writer.write("Allele,PeptideLength,ScoreBucket,Score,BucketCount,CumulativeCount");
            // writer.write(",CurrentScores,CurrentSize,CurrentBracket,CurrentSizeTotal");
            writer.newLine();

            // score each against each allele and build up a percentiles for each
            for(Map.Entry<String,Map<Integer,List<ScoreDistributionData>>> alleleEntry : mAlleleScoresMap.entrySet())
            {
                final String allele = alleleEntry.getKey();
                final Map<Integer,List<ScoreDistributionData>> pepLenDistMap = alleleEntry.getValue();

                for(List<ScoreDistributionData> pepLenScoreDist : pepLenDistMap.values())
                {
                    for(ScoreDistributionData scoreDist : pepLenScoreDist)
                    {
                        writer.write(String.format("%s,%d,%f,%.4f,%d,%d",
                                scoreDist.Allele, scoreDist.PeptideLength, scoreDist.ScoreBucket, scoreDist.Score,
                                scoreDist.BucketCount, scoreDist.CumulativeCount));

                        //writer.write(String.format(",%d,%f,%f,%f", currentScores, currentSize, currentBracket, currentSizeTotal));
                        writer.newLine();
                    }
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
