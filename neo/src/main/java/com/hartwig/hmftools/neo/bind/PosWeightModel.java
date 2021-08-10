package com.hartwig.hmftools.neo.bind;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.pow;
import static java.lang.Math.round;

import static com.hartwig.hmftools.neo.bind.BindConstants.ALLELE_POS_MAPPING_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACID_COUNT;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACID_C_FREQ_ADJUST;
import static com.hartwig.hmftools.neo.bind.BindConstants.ENTROPY_ADJUST;
import static com.hartwig.hmftools.neo.bind.BindConstants.ENTROPY_FACTOR;
import static com.hartwig.hmftools.neo.bind.BindConstants.MIN_OBSERVED_AA_POS_FREQ;
import static com.hartwig.hmftools.neo.bind.BindConstants.REF_PEPTIDE_LEFT_FIXED_POS;
import static com.hartwig.hmftools.neo.bind.BindConstants.REF_PEPTIDE_LENGTH;

import static org.apache.commons.math3.util.FastMath.log;

import java.io.BufferedWriter;
import java.util.List;
import java.util.Map;


import com.hartwig.hmftools.neo.utils.AminoAcidFrequency;

public class PosWeightModel
{
    private final CalcConstants mConstants;
    private final AminoAcidFrequency mAminoAcidFrequency;
    private final BlosumMapping mBlosumMapping;
    private final HlaSequences mHlaSequences;
    private final NoiseModel mNoiseModel;
    private final GlobalWeights mGlobalWeights;

    public PosWeightModel(final CalcConstants calcConstants, final HlaSequences hlaSequences)
    {
        mConstants = calcConstants;
        mAminoAcidFrequency = new AminoAcidFrequency();
        mBlosumMapping = new BlosumMapping();
        mHlaSequences = hlaSequences;
        mNoiseModel = new NoiseModel(mAminoAcidFrequency, calcConstants.NoiseProbability, calcConstants.NoiseWeight);
        mGlobalWeights = new GlobalWeights(calcConstants, mAminoAcidFrequency);
    }

    public boolean noiseEnabled() { return mNoiseModel.enabled(); }
    public final NoiseModel noiseModel() { return mNoiseModel; }

    public static final int INVALID_POS = -1;

    public final GlobalWeights getGlobalWeights() { return mGlobalWeights; }

    public static int peptidePositionToRef(int refLength, int peptideLength, int position)
    {
        if(position <= REF_PEPTIDE_LEFT_FIXED_POS)
            return position;

        // padded from the right
        return position + (refLength - peptideLength);
    }

    public static int refPeptidePositionToActual(int refLength, int peptideLength, int refPosition)
    {
        if(refPosition <= REF_PEPTIDE_LEFT_FIXED_POS)
            return refPosition;

        // padded from the right
        int actualPos = refPosition - (refLength - peptideLength);
        return actualPos > REF_PEPTIDE_LEFT_FIXED_POS ? actualPos : INVALID_POS;
    }

    public void buildNoiseCounts(final BindCountData bindCounts)
    {
        if(!mNoiseModel.enabled())
            return;

        final double[][] counts = bindCounts.getBindCounts();
        final double[][] noiseCounts = bindCounts.getNoiseCounts();

        for(int pos = 0; pos < bindCounts.PeptideLength; ++pos)
        {
            int posTotalBinds = 0;

            for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
            {
                posTotalBinds += counts[aa][pos];
            }

            if(posTotalBinds == 0)
                continue;

            for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
            {
                double obsCount = counts[aa][pos];
                double noiseCount = mNoiseModel.getExpected(aa, posTotalBinds, (int)round(obsCount));
                noiseCounts[aa][pos] = max(noiseCount, obsCount);
            }
        }
    }

    public void buildWeightedCounts(final BindCountData bindCounts, final List<BindCountData> peptideLengthCounts)
    {
        // translate the counts from various peptide lengths into this set of counts, normalising to the reference peptide length

        final double[][] weightedCounts = bindCounts.getWeightedCounts();

        for(BindCountData otherBindCounts : peptideLengthCounts)
        {
            final double[][] otherCounts = mNoiseModel.enabled() ? otherBindCounts.getNoiseCounts() : otherBindCounts.getBindCounts();
            int otherPeptideLength = otherBindCounts.PeptideLength;

            // use CountWeight for other peptide lengths' counts
            // CountWeight = 1 (PeptideLengthDiff) * 1 / (1 + (observedCount / LHW) ^ WeightExponent)
            double weight = 1;

            if(otherPeptideLength != bindCounts.PeptideLength)
            {
                weight = 1.0 / abs(bindCounts.PeptideLength - otherPeptideLength)
                        * 1 / (1 + pow(bindCounts.totalBindCount() / mConstants.PeptideLengthWeight, mConstants.WeightExponent));
            }

            for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
            {
                for(int otherPos = 0; otherPos < otherPeptideLength; ++otherPos)
                {
                    int refPos = peptidePositionToRef(REF_PEPTIDE_LENGTH, otherPeptideLength, otherPos);
                    int pos = refPeptidePositionToActual(REF_PEPTIDE_LENGTH, bindCounts.PeptideLength, refPos);

                    if(pos == INVALID_POS)
                        continue;

                    weightedCounts[aa][pos] += weight * otherCounts[aa][otherPos];
                }
            }
        }

        mGlobalWeights.processBindCounts(bindCounts);
    }

    public void buildFinalWeightedCounts(
            final BindCountData bindCounts, final List<BindCountData> allBindCounts, final Map<String,Integer> alleleTotalCounts)
    {
        // calculate Blosum similarity at this position vs all the other alleles from matching peptide lengths

        // WCount(A,L,P,AA) = LWCount(A,L,P,AA) + SUM(a<>A)  [ LWCount(a,L,P,AA)
        // * (2^(LogSim(m,M)) /  MAX(i=all motifs)[2^(LogSim(i,M))]] * [1 / ( 1+ Obs(A,L)/MHW)^E)]

        final double[][] weightedCounts = bindCounts.getWeightedCounts();
        final double[][] finalWeightedCounts = bindCounts.getFinalWeightedCounts();

        int alleleTotalCount = alleleTotalCounts.get(bindCounts.Allele);

        for(int pos = 0; pos < bindCounts.PeptideLength; ++pos)
        {
            // pos needs to be mapped to the 9-mer allele position peptide-position mapping
            int refPos = peptidePositionToRef(REF_PEPTIDE_LENGTH, bindCounts.PeptideLength, pos);
            int mappingPos = refPeptidePositionToActual(REF_PEPTIDE_LENGTH, ALLELE_POS_MAPPING_PEPTIDE_LENGTH, refPos);

            if(mappingPos == INVALID_POS)
            {
                // cannot use other alleles
                for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
                {
                    finalWeightedCounts[aa][pos] = weightedCounts[aa][pos];
                }

                continue;
            }

            String positionMotif = mHlaSequences.getSequence(bindCounts.Allele, mappingPos);

            if(positionMotif == null)
                continue;

            double selfScore = mBlosumMapping.calcSequenceBlosumScore(positionMotif);

            for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
            {
                for(BindCountData otherBindCounts : allBindCounts)
                {
                    if(otherBindCounts.Allele.equals(bindCounts.Allele))
                    {
                        // taken as is, no weighting
                        finalWeightedCounts[aa][pos] = weightedCounts[aa][pos];
                        continue;
                    }

                    double otherCount = otherBindCounts.getWeightedCounts()[aa][pos];

                    if(otherCount == 0)
                        continue;

                    String otherPositionMotif = mHlaSequences.getSequence(otherBindCounts.Allele, mappingPos);

                    if(otherPositionMotif == null)
                        continue;

                    double motifSimilarity = 1;

                    if(!otherPositionMotif.equals(positionMotif))
                    {
                        double crossAlleleScore = mBlosumMapping.calcSequenceBlosumScore(positionMotif, otherPositionMotif);
                        motifSimilarity = crossAlleleScore / selfScore;
                    }

                    double observationsWeight = 1 / (1 + pow(alleleTotalCount / mConstants.AlleleWeight, mConstants.WeightExponent));

                    double otherWeightedCount = otherCount * observationsWeight * motifSimilarity;
                    finalWeightedCounts[aa][pos] += otherWeightedCount;
                }
            }
        }
    }

    private static double[] calcPositionEntropy(final double[][] counts, final int peptideLength)
    {
        double[] entropy = new double[peptideLength];

        for(int pos = 0; pos < peptideLength; ++pos)
        {
            double posTotalCount = 0;

            for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
            {
                posTotalCount += counts[aa][pos];
            }

            double entropyTotal = 0;

            for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
            {
                double countsPerc = counts[aa][pos] / posTotalCount;

                if(countsPerc > 0)
                    entropyTotal -= countsPerc * log(2, countsPerc);
            }

            entropy[pos] = entropyTotal;
        }

        return entropy;
    }

    public BindScoreMatrix createMatrix(final BindCountData bindCounts)
    {
        mGlobalWeights.buildMatrixData();

        final double[][] finalWeightedCounts = bindCounts.getFinalWeightedCounts();

        BindScoreMatrix matrix = new BindScoreMatrix(bindCounts.Allele, bindCounts.PeptideLength);
        final double[][] data = matrix.getBindScores();

        /*
        double globalWeight = mConstants.GlobalWeight;

        double globalReductionFactor = 1;
        double[][] globalCounts = mGlobalWeights.get(bindCounts.PeptideLength);

        if(globalWeight > 0)
        {
            double total = calcCountsTotal(bindCounts.getFinalWeightedCounts());
            double globalTotal = mGlobalPeptideLengthTotals.get(bindCounts.PeptideLength);
            globalReductionFactor = total / globalTotal;
        }

        final double[] entropy = calcPositionEntropy(bindCounts.getWeightedCounts(), bindCounts.PeptideLength);
        */

        for(int pos = 0; pos < bindCounts.PeptideLength; ++pos)
        {
            double posTotalCount = 0;

            for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
            {
                posTotalCount += finalWeightedCounts[aa][pos];
            }

            // double entropyAdjust = ENTROPY_ADJUST > 0 ? 1 - pow(entropy[pos] / ENTROPY_FACTOR, ENTROPY_ADJUST) : 1;

            for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
            {
                double aaFrequency = mAminoAcidFrequency.getAminoAcidFrequency(aa);

                // Peptide Weight = log(max(posWeight(aa,pos) * AaAdjust, 0.005 * posWeightTotal)/AaFreq,2)
                double adjustedCount = finalWeightedCounts[aa][pos];

                if(aa == 'C')
                    adjustedCount *= AMINO_ACID_C_FREQ_ADJUST;

                // handle very low observation counts
                adjustedCount = max(adjustedCount, posTotalCount * MIN_OBSERVED_AA_POS_FREQ);

                /*
                if(globalWeight > 0)
                {
                    adjustedCount = globalWeight * globalReductionFactor * globalCounts[aa][pos] + (1 - globalWeight) * adjustedCount;
                }

                double posWeight = entropyAdjust * log(2, adjustedCount / aaFrequency);
                */

                double posWeight = log(2, adjustedCount / aaFrequency);
                data[aa][pos] = posWeight;
            }
        }

        return matrix;
    }

}
