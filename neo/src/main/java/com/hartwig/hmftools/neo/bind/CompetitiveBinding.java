package com.hartwig.hmftools.neo.bind;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Doubles;

public class CompetitiveBinding
{
    private final List<Double> mScoreRankBuckets;

    private final BufferedWriter mWriter;

    private static final double MIN_BUCKET_RANK = 0.00005;

    public CompetitiveBinding(final String outputFilename)
    {
        mScoreRankBuckets = Lists.newArrayListWithExpectedSize(16);

        double rankPercBucket = MIN_BUCKET_RANK;
        while(rankPercBucket < 1)
        {
            mScoreRankBuckets.add(rankPercBucket);
            rankPercBucket *= 2;
        }

        mWriter = initWriter(outputFilename);
    }

    public void processAllelePeptideRanks(final Map<String,Map<Integer,List<BindData>>> allelePeptideData)
    {
        for(Map.Entry<String, Map<Integer, List<BindData>>> alleleEntry : allelePeptideData.entrySet())
        {
            final String allele = alleleEntry.getKey();
            final Map<Integer,List<BindData>> pepLenBindDataMap = alleleEntry.getValue();

            final Map<Integer,List<Double>> pepLenRanks = Maps.newHashMap();

            for(Map.Entry<Integer,List<BindData>> pepLenEntry : pepLenBindDataMap.entrySet())
            {
                int peptideLength = pepLenEntry.getKey();
                List<Double> peptideRanks = Lists.newArrayListWithCapacity(pepLenEntry.getValue().size());
                pepLenEntry.getValue().forEach(x -> peptideRanks.add(x.rankPercentile()));
                pepLenRanks.put(peptideLength, peptideRanks);
            }

            processAllelePeptideRanks(allele, pepLenRanks);
        }

        closeBufferedWriter(mWriter);
    }

    private void processAllelePeptideRanks(final String allele, final Map<Integer,List<Double>> pepLenRanks)
    {
        Map<Integer,int[]> pepLenRankCounts = Maps.newHashMap();
        int totalBuckets = mScoreRankBuckets.size();

        for(Map.Entry<Integer,List<Double>> entry : pepLenRanks.entrySet())
        {
            int peptideLength = entry.getKey();

            int[] rankCounts = new int[totalBuckets];
            pepLenRankCounts.put(peptideLength, rankCounts);

            for(double rank : entry.getValue())
            {
                for(int i = 0; i < totalBuckets; ++i)
                {
                    double bucketRank = mScoreRankBuckets.get(i);

                    if(rank > bucketRank)
                        continue;

                    if(rank < bucketRank || Doubles.equal(bucketRank, rank))
                    {
                        ++rankCounts[i];
                        break;
                    }

                    if(i < totalBuckets - 1)
                    {
                        double nextBucketRank = mScoreRankBuckets.get(i + 1);

                        if(rank < nextBucketRank || Doubles.equal(nextBucketRank, rank))
                        {
                            ++rankCounts[i + 1];
                            break;
                        }
                    }
                    else
                    {
                        ++rankCounts[totalBuckets - 1];
                    }
                }
            }
        }

        Map<Integer,double[]> pepLenLikelihoods = Maps.newHashMap();
        double totalAdjustedCounts = 0;

        for(Map.Entry<Integer,List<Double>> entry : pepLenRanks.entrySet())
        {
            int peptideLength = entry.getKey();

            double[] likelihoods = new double[totalBuckets];
            pepLenLikelihoods.put(peptideLength, likelihoods);

            final int[] rankCounts = pepLenRankCounts.get(peptideLength);

            for(int i = rankCounts.length - 1; i >= 0; --i)
            {
                // MAX(D4/($B4/0.00005),IF($B5>$B4,I5,0))
                double adjCount = rankCounts[i] / (mScoreRankBuckets.get(i) / MIN_BUCKET_RANK);

                if(i == rankCounts.length - 1)
                {
                    likelihoods[i] = adjCount;
                }
                else
                {
                    // will always be at least as high as the next bucket up
                    likelihoods[i] = max(adjCount, likelihoods[i + 1]);
                }

                totalAdjustedCounts += likelihoods[i];
            }
        }

        for(Map.Entry<Integer,List<Double>> entry : pepLenRanks.entrySet())
        {
            int peptideLength = entry.getKey();

            double[] likelihoods = pepLenLikelihoods.get(peptideLength);

            for(int i = 0; i < likelihoods.length; ++i)
            {
                likelihoods[i] /= totalAdjustedCounts;
            }
        }

        writeAlleleData(allele, pepLenLikelihoods);
    }

    private BufferedWriter initWriter(final String filename)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Allele,PeptideLength");

            for(Double rankBucket : mScoreRankBuckets)
            {
                writer.write(String.format(",P_%f", rankBucket));
            }

            writer.newLine();

            return writer;
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to initialise relative presentation file file({}): {}", filename, e.toString());
            return null;
        }
    }

    private void writeAlleleData(final String allele, final Map<Integer,double[]> pepLenLikelihoodMap)
    {
        if(mWriter == null)
            return;

        try
        {
            for(Map.Entry<Integer,double[]> entry : pepLenLikelihoodMap.entrySet())
            {
                int peptideLength = entry.getKey();
                final double[] likelihoods = entry.getValue();

                mWriter.write(String.format("%s,%d", allele, peptideLength));

                for(int i = 0; i < likelihoods.length; ++i)
                {
                    mWriter.write(String.format(",%.6f", likelihoods[i]));
                }

                mWriter.newLine();
            }
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to initialise relative presentation file: {}", e.toString());
        }
    }



}
