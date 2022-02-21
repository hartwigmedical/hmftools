package com.hartwig.hmftools.neo.bind;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.DELIM;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_ALLELE;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PEPTIDE_LEN;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_PEPTIDE_LENGTHS;
import static com.hartwig.hmftools.neo.bind.BindConstants.MIN_LIKELIHOOD_ALLELE_BIND_COUNT;
import static com.hartwig.hmftools.neo.bind.BindConstants.MIN_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.BindConstants.REF_PEPTIDE_LENGTH;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Doubles;

public class BindingLikelihood
{
    // statically-defined rank buckets in exponentially increasing size
    private final List<Double> mScoreRankBuckets;

    // map of allele to a matrix of peptide-length and rank buckets, containing the likelihood values
    private final Map<String,double[][]> mAlleleLikelihoodMap;

    private BufferedWriter mWriter;
    private boolean mHasData;

    private static final double INVALID_LIKELIHOOD = -1;
    private static final double MIN_BUCKET_RANK = 0.00005;
    private static final int PEPTIDE_LENGTHS = 5;

    private static final double MIN_EMPTY_PEP_LEN_LIKELIHOOD = 0.25;
    private static final double MIN_EMPTY_PEP_LEN_FACTOR = 1000;
    private static final double MIN_EMPTY_PEP_LEN_BUCKET_THRESHOLD = 0.01;

    public BindingLikelihood()
    {
        mWriter = null;
        mAlleleLikelihoodMap = Maps.newHashMap();
        mHasData = false;

        mScoreRankBuckets = Lists.newArrayListWithExpectedSize(16);

        double rankPercBucket = MIN_BUCKET_RANK;
        while(rankPercBucket < 1)
        {
            mScoreRankBuckets.add(rankPercBucket);
            rankPercBucket *= 2;
        }
    }

    public boolean hasData() { return mHasData; }

    public double getBindingLikelihood(final String allele, final String peptide, final double rank)
    {
        if(!mHasData)
            return INVALID_LIKELIHOOD;

        int peptideLength = peptide.length();
        int pepLenIndex = peptideLengthIndex(peptideLength);

        if(pepLenIndex == INVALID_PEP_LEN)
            return INVALID_LIKELIHOOD;

        double[][] likelihoods = mAlleleLikelihoodMap.get(allele);
        if(likelihoods == null)
        {
            likelihoods = mAlleleLikelihoodMap.get(extractTwoDigitType(allele));

            if(likelihoods == null)
                likelihoods = mAlleleLikelihoodMap.get(extractGene(allele));

            if(likelihoods == null)
                return INVALID_LIKELIHOOD;
        }

        for(int i = 0; i < mScoreRankBuckets.size(); ++i)
        {
            if(rank < mScoreRankBuckets.get(i))
            {
                double likelihood = likelihoods[pepLenIndex][i];

                if(i == 0)
                    return likelihood;

                double lowerLikelihood = likelihoods[pepLenIndex][i - 1];
                double lowerRank = mScoreRankBuckets.get(i - 1);
                double upperRank = mScoreRankBuckets.get(i);
                double upperPerc = (rank - lowerRank) / (upperRank - lowerRank);
                return upperPerc * likelihood + (1 - upperPerc) * lowerLikelihood;
            }
        }

        return 0;
    }

    public boolean loadLikelihoods(final String filename)
    {
        if(filename == null)
            return false;

        try
        {
            final List<String> lines = Files.readAllLines(new File(filename).toPath());

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIM);
            lines.remove(0);

            int alleleIndex = fieldsIndexMap.get(FLD_ALLELE);
            int pepLenIndex = fieldsIndexMap.get(FLD_PEPTIDE_LEN);

            String currentAllele = "";
            double[][] likelihoods = null;

            for(String line : lines)
            {
                String[] items = line.split(DELIM, -1);
                String allele = items[alleleIndex];
                int peptideLength = Integer.parseInt(items[pepLenIndex]);

                if(!currentAllele.equals(allele))
                {
                    currentAllele = allele;
                    likelihoods = new double[PEPTIDE_LENGTHS][mScoreRankBuckets.size()];
                    mAlleleLikelihoodMap.put(allele, likelihoods);
                }

                int index = pepLenIndex + 1;
                for(int i = 0; index < items.length; ++i, ++index)
                {
                    likelihoods[peptideLengthIndex(peptideLength)][i] = Double.parseDouble(items[index]);
                }
            }

            NE_LOGGER.info("loaded {} alleles peptide likelihoods from {}", mAlleleLikelihoodMap.size(), filename);
            mHasData = true;
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to read peptide likelihoods file: {}", e.toString());
            return false;
        }

        return true;
    }

    private static int INVALID_PEP_LEN = -1;

    private static int peptideLengthIndex(int peptideLength)
    {
        if(peptideLength < MIN_PEPTIDE_LENGTH || peptideLength > REF_PEPTIDE_LENGTH)
            return INVALID_PEP_LEN;

        return peptideLength - MIN_PEPTIDE_LENGTH;
    }

    private static int peptideLengthFromIndex(int index)
    {
        return MIN_PEPTIDE_LENGTH + index;
    }

    public void buildAllelePeptideLikelihoods(
            final Map<String,Map<Integer,List<BindData>>> allelePeptideData, final String outputFilename)
    {
        if(outputFilename != null)
            mWriter = initWriter(outputFilename);

        // any allele with sufficient counts will have a likelihood distribution calculated for it

        // in addition, distributions will be calculated for 2-digit alleles and at the HLA gene level
        final Map<String,Map<Integer,List<Double>>> sharedPepLenRanks = Maps.newHashMap();

        for(Map.Entry<String, Map<Integer, List<BindData>>> alleleEntry : allelePeptideData.entrySet())
        {
            final String allele = alleleEntry.getKey();
            final Map<Integer,List<BindData>> pepLenBindDataMap = alleleEntry.getValue();

            int alleleBindCount = pepLenBindDataMap.values().stream().mapToInt(x -> x.size()).sum();

            if(alleleBindCount < MIN_LIKELIHOOD_ALLELE_BIND_COUNT)
                continue;

            String hlaGene = extractGene(allele);
            String twoDigitType = extractTwoDigitType(allele);

            Map<Integer,List<Double>> genePepLenRanks = sharedPepLenRanks.get(hlaGene);
            if(genePepLenRanks == null)
            {
                genePepLenRanks = Maps.newHashMap();
                sharedPepLenRanks.put(hlaGene, genePepLenRanks);
            }

            Map<Integer,List<Double>> twoDigitPepLenRanks = sharedPepLenRanks.get(twoDigitType);
            if(twoDigitPepLenRanks == null)
            {
                twoDigitPepLenRanks = Maps.newHashMap();
                sharedPepLenRanks.put(twoDigitType, twoDigitPepLenRanks);
            }

            final Map<Integer,List<Double>> pepLenRanks = Maps.newHashMap();

            for(Map.Entry<Integer,List<BindData>> pepLenEntry : pepLenBindDataMap.entrySet())
            {
                int peptideLength = pepLenEntry.getKey();
                List<Double> peptideRanks = Lists.newArrayListWithCapacity(pepLenEntry.getValue().size());
                pepLenEntry.getValue().forEach(x -> peptideRanks.add(x.rankPercentile()));
                pepLenRanks.put(peptideLength, peptideRanks);

                List<Double> sharedRanks = genePepLenRanks.get(peptideLength);
                if(sharedRanks == null)
                {
                    sharedRanks = Lists.newArrayList();
                    genePepLenRanks.put(peptideLength, sharedRanks);
                }

                sharedRanks.addAll(peptideRanks);

                sharedRanks = twoDigitPepLenRanks.get(peptideLength);
                if(sharedRanks == null)
                {
                    sharedRanks = Lists.newArrayList();
                    twoDigitPepLenRanks.put(peptideLength, sharedRanks);
                }

                sharedRanks.addAll(peptideRanks);
            }

            buildAllelePeptideLikelihoods(allele, pepLenRanks);
        }

        for(Map.Entry<String,Map<Integer,List<Double>>> sharedAlleleEntry : sharedPepLenRanks.entrySet())
        {
            buildAllelePeptideLikelihoods(sharedAlleleEntry.getKey(), sharedAlleleEntry.getValue());
        }

        mHasData = true;

        closeBufferedWriter(mWriter);
    }

    private static String extractGene(final String allele)
    {
        // A, B or C
        return allele.substring(0, 1);
    }

    private static String extractTwoDigitType(final String allele)
    {
        // eg A02
        return allele.substring(0, 3);
    }

    private void buildAllelePeptideLikelihoods(final String allele, final Map<Integer,List<Double>> pepLenRanks)
    {
        Map<Integer,double[]> pepLenRankCounts = Maps.newHashMap();
        int totalBuckets = mScoreRankBuckets.size();
        int totalPositivesCount = 0;

        for(int peptideLength : DEFAULT_PEPTIDE_LENGTHS)
        {
            double[] rankCounts = new double[totalBuckets];
            pepLenRankCounts.put(peptideLength, rankCounts);

            List<Double> pepLengthRanks = pepLenRanks.get(peptideLength);

            if(pepLengthRanks == null)
                continue;

            totalPositivesCount += pepLengthRanks.size();

            for(double rank : pepLengthRanks)
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

        // fill in values for zeros using a fraction of total positives for lengths with none, and halving for other missing buckets
        double minLikelihood = min(MIN_EMPTY_PEP_LEN_LIKELIHOOD, totalPositivesCount / MIN_EMPTY_PEP_LEN_FACTOR);

        for(int peptideLength : DEFAULT_PEPTIDE_LENGTHS)
        {
            double[] rankCounts = pepLenRankCounts.get(peptideLength);

            for(int i = 0; i < rankCounts.length; ++i)
            {
                if(rankCounts[i] > 0)
                    continue;

                double bucketRank = mScoreRankBuckets.get(i);
                if(bucketRank < MIN_EMPTY_PEP_LEN_BUCKET_THRESHOLD)
                {
                    rankCounts[i] = minLikelihood;
                }
                else
                {
                    rankCounts[i] = min(minLikelihood, rankCounts[i - 1] * 0.5);
                }
            }
        }

        double totalAdjustedCounts = 0;

        double[][] pepLenLikelihoods = new double[PEPTIDE_LENGTHS][totalBuckets];

        for(int peptideLength : DEFAULT_PEPTIDE_LENGTHS)
        {
            int pepLenIndex = peptideLengthIndex(peptideLength);

            final double[] rankCounts = pepLenRankCounts.get(peptideLength);

            for(int i = rankCounts.length - 1; i >= 0; --i)
            {
                double adjCount = rankCounts[i] / (mScoreRankBuckets.get(i) / MIN_BUCKET_RANK);

                if(i == rankCounts.length - 1)
                {
                    pepLenLikelihoods[pepLenIndex][i] = adjCount;
                }
                else
                {
                    // will always be at least as high as the next bucket up
                    pepLenLikelihoods[pepLenIndex][i] = max(adjCount, pepLenLikelihoods[pepLenIndex][i + 1]);
                }

                totalAdjustedCounts += pepLenLikelihoods[pepLenIndex][i];
            }
        }

        for(int pl = 0; pl < pepLenLikelihoods.length; ++pl)
        {
            for(int i = 0; i < totalBuckets; ++i)
            {
                pepLenLikelihoods[pl][i] /= totalAdjustedCounts;
            }
        }

        mAlleleLikelihoodMap.put(allele, pepLenLikelihoods);

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
            NE_LOGGER.error("failed to initialise binding likelihood file file({}): {}", filename, e.toString());
            return null;
        }
    }

    private void writeAlleleData(final String allele, final double[][] pepLenLikelihoods)
    {
        if(mWriter == null)
            return;

        try
        {
            for(int pl = 0; pl < pepLenLikelihoods.length; ++pl)
            {
                mWriter.write(String.format("%s,%d", allele, peptideLengthFromIndex(pl)));

                for(int i = 0; i < mScoreRankBuckets.size(); ++i)
                {
                    mWriter.write(String.format(",%4.3e", pepLenLikelihoods[pl][i]));
                }

                mWriter.newLine();
            }
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to initialise binding likelihood file: {}", e.toString());
        }
    }



}
