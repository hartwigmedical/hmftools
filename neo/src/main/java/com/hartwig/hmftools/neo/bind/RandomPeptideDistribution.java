package com.hartwig.hmftools.neo.bind;

import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;

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

    private final Map<String,List<ScoreDistributionData>> mAlleleScoresMap;

    public RandomPeptideDistribution(final BinderConfig config)
    {
        mInputPeptidesFilename = config.RandomPeptidesFile;
        mDistributionFilename = config.RandomPeptideDistributionFile;

        mAlleleScoresMap = Maps.newHashMap();

        mDataLoaded = loadDistribution();
    }

    public boolean hasData() { return mDataLoaded; }

    public Map<String,List<ScoreDistributionData>> getAlleleScoresMap() { return mAlleleScoresMap; }

    public double getScoreRank(final String allele, double score)
    {
        List<ScoreDistributionData> scores = mAlleleScoresMap.get(allele);

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
        if(!Files.exists(Paths.get(mDistributionFilename)))
            return false;

        try
        {
            final List<String> lines = Files.readAllLines(new File(mDistributionFilename).toPath());

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            for(String line : lines)
            {
                ScoreDistributionData data = ScoreDistributionData.fromCsv(line, fieldsIndexMap);

                List<ScoreDistributionData> dataList = mAlleleScoresMap.get(data.Allele);

                if(dataList == null)
                {
                    dataList = Lists.newArrayList();
                    mAlleleScoresMap.put(data.Allele, dataList);
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

    public void buildDistribution(final List<BindScoreMatrix> matrixList)
    {
        List<String> randomPeptides = Lists.newArrayList();

        try
        {
            randomPeptides.addAll(Files.readAllLines(new File(mInputPeptidesFilename).toPath()));
            randomPeptides.remove(0);

            NE_LOGGER.info("loaded {} random peptides", randomPeptides.size());

            BufferedWriter writer = createBufferedWriter(mDistributionFilename, false);

            writer.write("Allele,PeptideLength,ScoreBucket,Score");
            // writer.write(",CurrentScores,CurrentSize,CurrentBracket,CurrentSizeTotal");
            writer.newLine();

            int peptideCount = randomPeptides.size();

            // score each against each allele and build up a percentiles for each
            for(BindScoreMatrix matrix : matrixList)
            {
                NE_LOGGER.info("building distribution for allele({}) peptideLength({})",
                        matrix.Allele, matrix.PeptideCount);

                List<Double> peptideScores = Lists.newArrayList();

                int count = 0;
                for(String peptide : randomPeptides)
                {
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
                int requiredScores = (int)round(peptideCount * currentSize);

                double scoreTotal = 0;
                int currentScores = 0;
                double currentSizeTotal = 0;

                for(Double score : peptideScores)
                {
                    scoreTotal += score;
                    ++currentScores;

                    if(currentScores >= requiredScores)
                    {
                        double avgScore = scoreTotal / currentScores;

                        writer.write(String.format("%s,%d,%f,%.4f", matrix.Allele, matrix.PeptideCount, currentSizeTotal, avgScore));

                        //writer.write(String.format(",%d,%f,%f,%f", currentScores, currentSize, currentBracket, currentSizeTotal));
                        writer.newLine();

                        scoreTotal = 0;
                        currentScores = 0;
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

            writer.close();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write random peptide file: {}", e.toString());
        }
    }

}
