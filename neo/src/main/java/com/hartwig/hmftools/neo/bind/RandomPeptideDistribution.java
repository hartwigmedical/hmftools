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
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sigs.VectorUtils;
import com.hartwig.hmftools.common.utils.Doubles;

public class RandomPeptideDistribution
{
    private final String mInputPeptidesFilename;
    private final String mDistributionFilename;
    private final boolean mDataLoaded;

    private final Map<String,Map<Integer,List<ScoreDistributionData>>> mAlleleScoresMap; // allele to peptide length to distribution

    public RandomPeptideDistribution(final BinderConfig config)
    {
        mInputPeptidesFilename = config.RandomPeptidesFile;

        mAlleleScoresMap = Maps.newHashMap();

        if(config.RandomPeptideDistributionFile != null)
        {
            mDistributionFilename = config.RandomPeptideDistributionFile;
            mDataLoaded = loadDistribution();
        }
        else
        {
            mDistributionFilename = config.formFilename("random_peptide_dist");
            mDataLoaded = false;
        }
    }

    public boolean hasData() { return mDataLoaded; }

    public Map<String,Map<Integer,List<ScoreDistributionData>>> getAlleleScoresMap() { return mAlleleScoresMap; }

    public double getScoreRank(final String allele, final int peptideLength, double score)
    {
        Map<Integer,List<ScoreDistributionData>> peptideLengthMap = mAlleleScoresMap.get(allele);

        if(peptideLengthMap == null)
            return -1;

        List<ScoreDistributionData> scores = peptideLengthMap.get(peptideLength);

        if(scores == null || scores.isEmpty())
            return -1;

        if(score > scores.get(0).Score)
            return scores.get(0).ScoreBucket * 0.5;

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

    private boolean loadDistribution()
    {
        if(mDistributionFilename == null || !Files.exists(Paths.get(mDistributionFilename)))
            return false;

        try
        {
            final List<String> lines = Files.readAllLines(new File(mDistributionFilename).toPath());

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
            NE_LOGGER.error("failed to read training binding data file: {}", e.toString());
            return false;
        }

        return true;
    }

    public void buildDistribution(final Map<String,Map<Integer,BindScoreMatrix>> alleleBindMatrixMap)
    {
        List<String> refRandomPeptides = Lists.newArrayList();

        if(mInputPeptidesFilename != null)
        {
            try
            {
                refRandomPeptides.addAll(Files.readAllLines(new File(mInputPeptidesFilename).toPath()));
                refRandomPeptides.remove(0);

                NE_LOGGER.info("loaded {} random peptides", refRandomPeptides.size());
            }
            catch(IOException e)
            {
                NE_LOGGER.error("failed to write random peptide file: {}", e.toString());
                return;
            }
        }
        else
        {
            // take the random peptides from any of the provided alleles??
            /*
            randomPeptides = alleleBindDataMap.get(matrix.Allele).stream()
                .filter(x -> x.isRandom())
                .filter(x -> x.peptideLength() == matrix.PeptideLength)
                .map(x -> x.Peptide).collect(Collectors.toList());
             */
        }

        try
        {
            BufferedWriter writer = createBufferedWriter(mDistributionFilename, false);

            writer.write("Allele,PeptideLength,ScoreBucket,Score,BucketCount,CumulativeCount");
            // writer.write(",CurrentScores,CurrentSize,CurrentBracket,CurrentSizeTotal");
            writer.newLine();

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

                    List<String> randomPeptides = Lists.newArrayList(refRandomPeptides);

                    int peptideCount = randomPeptides.size();

                    List<Double> peptideScores = Lists.newArrayList();

                    int count = 0;

                    for(String peptide : randomPeptides)
                    {
                        if(peptide.length() > matrix.PeptideLength)
                            peptide = peptide.substring(0, matrix.PeptideLength);

                        double score = matrix.calcScore(peptide);
                        VectorUtils.optimisedAdd(peptideScores, score, false);

                        ++count;

                        if(count > 0 && (count % 100000) == 0)
                        {
                            NE_LOGGER.debug("added {} sorted random peptide scores", count);
                        }
                    }

                    // write the distribution as 0.0001 up to 0.01, 0.001 up to 0.01, then 0.01 up to 100%
                    List<Double> scoreDiscreteSizes = Lists.newArrayList(0.0001, 0.001, 0.01);
                    List<Double> scoreDiscreteBrackets = Lists.newArrayList(0.01, 0.1, 1.0);
                    int discreteIndex = 0;
                    double currentSize = scoreDiscreteSizes.get(discreteIndex);
                    double currentBracket = scoreDiscreteBrackets.get(discreteIndex);
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

                            writer.write(String.format("%s,%d,%f,%.4f,%d,%d",
                                    matrix.Allele, matrix.PeptideLength, currentSizeTotal, avgScore, currentScoreCount, cumulativeScores));

                            //writer.write(String.format(",%d,%f,%f,%f", currentScores, currentSize, currentBracket, currentSizeTotal));
                            writer.newLine();

                            scoresDistributions.add(new ScoreDistributionData(matrix.Allele, matrix.PeptideLength, currentSizeTotal, avgScore));

                            scoreTotal = 0;
                            currentScoreCount = 0;
                            currentSizeTotal += currentSize;

                            if(Doubles.equal(currentSizeTotal, currentBracket))
                            {
                                ++discreteIndex;

                                if(discreteIndex < scoreDiscreteSizes.size())
                                {
                                    currentSize = scoreDiscreteSizes.get(discreteIndex);
                                    currentBracket = scoreDiscreteBrackets.get(discreteIndex);
                                    requiredScores = (int) round(peptideCount * currentSize);
                                }
                            }
                        }
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

}
