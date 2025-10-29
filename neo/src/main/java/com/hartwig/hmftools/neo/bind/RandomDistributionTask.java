package com.hartwig.hmftools.neo.bind;

import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.PAN_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.BindConstants.PAN_PEPTIDE_MAX_LENGTH;
import static com.hartwig.hmftools.neo.bind.RandomPeptideDistribution.getScoreRank;

import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.VectorUtils;
import com.hartwig.hmftools.neo.PeptideData;

public class RandomDistributionTask implements Callable<Void>
{
    private final String mAllele;

    private final Map<Integer, List<PeptideData>> mRandomPeptideMap; // by length and with flanking data
    private final Map<String, Map<Integer, List<ScoreDistributionData>>> mAlleleScoreDistributions;
    private final FlankScores mFlankScores;
    private final Map<Integer, BindScoreMatrix> mPeptideLengthMatrixMap;

    private final int mTaskType;

    // outputs
    private final Map<Integer, List<ScoreDistributionData>> mPeptideLengthDistributions; // peptide length to distribution
    private final List<ScoreDistributionData> mLikelihoodDistributions; // distribution of likelihoods

    private final BindingLikelihood mBindingLikelihood;
    private final ExpressionLikelihood mExpressionLikelihood;

    private final List<double[]> mDiscreteScoreData;

    private static final int SCORE_SIZE = 0;
    private static final int SCORE_BRACKET = 1;

    public static final int TASK_TYPE_SCORE_RANK = 0;
    public static final int TASK_TYPE_LIKELIHOOD_RANK = 1;

    // instantiate one of 2 tasks
    public RandomDistributionTask(
            final String allele, final Map<Integer, BindScoreMatrix> peptideLengthMatrixMap,
            final Map<Integer, List<PeptideData>> randomPeptideMap, final FlankScores flankScores)
    {
        this(TASK_TYPE_SCORE_RANK, allele, peptideLengthMatrixMap, randomPeptideMap, flankScores,
                null, null, null);
    }

    public RandomDistributionTask(
            final String allele, final Map<Integer, BindScoreMatrix> peptideLengthMatrixMap,
            final Map<Integer, List<PeptideData>> randomPeptideMap, final FlankScores flankScores,
            final Map<String, Map<Integer, List<ScoreDistributionData>>> alleleScoreDistributions,
            final BindingLikelihood bindingLikelihood, final ExpressionLikelihood expressionLikelihood)
    {
        this(TASK_TYPE_LIKELIHOOD_RANK, allele, peptideLengthMatrixMap, randomPeptideMap, flankScores,
                alleleScoreDistributions, bindingLikelihood, expressionLikelihood);
    }

    private RandomDistributionTask(
            final int taskType, final String allele, final Map<Integer, BindScoreMatrix> peptideLengthMatrixMap,
            final Map<Integer, List<PeptideData>> randomPeptideMap, final FlankScores flankScores,
            final Map<String, Map<Integer, List<ScoreDistributionData>>> alleleScoreDistributions,
            final BindingLikelihood bindingLikelihood, final ExpressionLikelihood expressionLikelihood)
    {
        mTaskType = taskType;
        mAllele = allele;
        mPeptideLengthMatrixMap = peptideLengthMatrixMap;
        mRandomPeptideMap = randomPeptideMap;
        mFlankScores = flankScores;
        mAlleleScoreDistributions = alleleScoreDistributions;
        mBindingLikelihood = bindingLikelihood;
        mExpressionLikelihood = expressionLikelihood;

        mPeptideLengthDistributions = Maps.newHashMap();
        mLikelihoodDistributions = Lists.newArrayList();

        mDiscreteScoreData = generateDistributionBuckets();
    }

    public static List<double[]> generateDistributionBuckets()
    {
        List<double[]> discreteScoreData = Lists.newArrayList();
        discreteScoreData.add(new double[] { 0.0001, 0.01 });
        discreteScoreData.add(new double[] { 0.001, 0.1 });
        discreteScoreData.add(new double[] { 0.01, 0.25 });
        discreteScoreData.add(new double[] { 0.05, 1.0 });
        return discreteScoreData;
    }

    public final String allele() { return mAllele; }

    public Map<Integer, List<ScoreDistributionData>> getPeptideLengthScoreDistributions() { return mPeptideLengthDistributions; }
    public List<ScoreDistributionData> getLikelihoodDistributions() { return mLikelihoodDistributions; }

    @Override
    public Void call()
    {
        if(mTaskType == TASK_TYPE_SCORE_RANK)
            buildScoreDistribution();
        else
            buildLikelihoodDistribution();

        return null;
    }

    private void buildScoreDistribution()
    {
        if(mRandomPeptideMap.isEmpty())
            return;

        // score each against each allele and build up a percentiles for each
        NE_LOGGER.debug("building distribution for allele({})", mAllele);

        for(BindScoreMatrix matrix : mPeptideLengthMatrixMap.values())
        {
            List<PeptideData> randomPeptides = mRandomPeptideMap.get(matrix.PeptideLength);

            if(randomPeptides == null || randomPeptides.isEmpty())
            {
                NE_LOGGER.error("missing random peptide data for length({})", matrix.PeptideLength);
                return;
            }

            List<Double> peptideScores = Lists.newArrayListWithExpectedSize(randomPeptides.size());

            int count = 0;

            for(PeptideData peptideData : randomPeptides)
            {
                double score = BindScorer.calcScore(matrix, mFlankScores, peptideData.Peptide, peptideData.UpFlank, peptideData.DownFlank);

                VectorUtils.optimisedAdd(peptideScores, score, false);

                ++count;

                if(count > 0 && (count % 500000) == 0)
                {
                    NE_LOGGER.debug("added {} sorted random peptide scores", count);
                }
            }

            List<ScoreDistributionData> scoresDistributions = generateDistribution(matrix.Allele, matrix.PeptideLength, peptideScores);
            mPeptideLengthDistributions.put(matrix.PeptideLength, scoresDistributions);
        }
    }

    public void buildLikelihoodDistribution()
    {
        if(mRandomPeptideMap.isEmpty())
            return;

        NE_LOGGER.debug("building likelihood distribution for allele({})", mAllele);

        List<Double> likelihoodScores = Lists.newArrayList();

        for(BindScoreMatrix matrix : mPeptideLengthMatrixMap.values())
        {
            // for now hard-code exclusion of length 12 since other tools are 8-11
            if(matrix.PeptideLength > PAN_PEPTIDE_MAX_LENGTH)
                continue;

            int count = 0;

            List<PeptideData> randomPeptides = mRandomPeptideMap.get(matrix.PeptideLength);

            if(randomPeptides == null || randomPeptides.isEmpty())
                return;

            for(PeptideData peptideData : randomPeptides)
            {
                double score = BindScorer.calcScore(matrix, mFlankScores, peptideData.Peptide, peptideData.UpFlank, peptideData.DownFlank);
                double rank = getScoreRank(mAlleleScoreDistributions, mAllele, matrix.PeptideLength, score);
                double likelihood = mBindingLikelihood.getBindingLikelihood(mAllele, peptideData.Peptide, rank);

                if(likelihood > 0 && mExpressionLikelihood != null && mExpressionLikelihood.hasData())
                {
                    double expLikelihood = mExpressionLikelihood.calcLikelihood(peptideData.TPM);
                    likelihood *= expLikelihood;
                }

                VectorUtils.optimisedAdd(likelihoodScores, likelihood, false);

                ++count;

                if(count > 0 && (count % 500000) == 0)
                {
                    NE_LOGGER.debug("added {} sorted random peptide likelihood scores", count);
                }
            }
        }

        mLikelihoodDistributions.addAll(generateDistribution(mAllele, PAN_PEPTIDE_LENGTH, likelihoodScores));
    }

    private List<ScoreDistributionData> generateDistribution(final String allele, final int peptideLength, final List<Double> peptideScores)
    {
        return generateDistribution(allele, peptideLength, peptideScores, mDiscreteScoreData);
    }

    public static List<ScoreDistributionData> generateDistribution(
            final String allele, final int peptideLength, final List<Double> peptideScores, final List<double[]> discreteScoreData)
    {
        if(peptideScores.size() < 2)
            return Lists.newArrayList();

        // write the distribution frequency by grouping the sorted items into increasing discrete bucket increments
        // store the best score at each point
        int totalScores = peptideScores.size();
        int discreteIndex = 0;
        double currentSize = discreteScoreData.get(discreteIndex)[SCORE_SIZE];
        double currentBracket = discreteScoreData.get(discreteIndex)[SCORE_BRACKET];

        int requiredScores = (int) max(round(totalScores * currentSize), 1);

        int currentScoreCount = 0;
        double scoreTotal = 0;

        List<ScoreDistributionData> scoresDistributions = Lists.newArrayList();

        // start with the top score in the first bucket by itself
        scoresDistributions.add(new ScoreDistributionData(
                allele, peptideLength, scoreTotal, peptideScores.get(0), 1, 1));

        currentScoreCount = 0;
        scoreTotal += currentSize;

        for(int i = 1; i < peptideScores.size() - 1; ++i)
        {
            Double score = peptideScores.get(i);
            ++currentScoreCount;

            if(currentScoreCount >= requiredScores)
            {
                scoresDistributions.add(new ScoreDistributionData(allele, peptideLength, scoreTotal, score, currentScoreCount, i + 1));

                currentScoreCount = 0;

                if(Doubles.equal(scoreTotal, currentBracket))
                {
                    ++discreteIndex;

                    if(discreteIndex < discreteScoreData.size())
                    {
                        currentSize = discreteScoreData.get(discreteIndex)[SCORE_SIZE];
                        currentBracket = discreteScoreData.get(discreteIndex)[SCORE_BRACKET];
                        requiredScores = (int) max(round(totalScores * currentSize), 1);
                    }
                }

                scoreTotal += currentSize;
            }
        }

        // add the last score so as to interpolate between scores right up to 100%
        double lastScore = peptideScores.get(peptideScores.size() - 1);
        scoresDistributions.add(new ScoreDistributionData(
                allele, peptideLength, 1, lastScore, 1, peptideScores.size()));

        return scoresDistributions;
    }
}
