package com.hartwig.hmftools.neo.bind;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.pow;
import static java.lang.Math.round;

import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.ALLELE_POS_MAPPING_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACIDS;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACID_COUNT;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACID_C_FREQ_ADJUST;
import static com.hartwig.hmftools.neo.bind.BindConstants.MIN_OBSERVED_AA_POS_FREQ;
import static com.hartwig.hmftools.neo.bind.BindConstants.REF_PEPTIDE_LEFT_FIXED_POS;
import static com.hartwig.hmftools.neo.bind.BindConstants.REF_PEPTIDE_LENGTH;

import static org.apache.commons.math3.util.FastMath.log;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.neo.utils.AminoAcidFrequency;

public class PosWeightModel
{
    private final CalcConstants mConstants;
    private final AminoAcidFrequency mAminoAcidFrequency;
    private final BlosumMapping mBlosumMapping;
    private final HlaSequences mHlaSequences;
    private final NoiseModel mNoiseModel;

    private final Map<Integer,double[][]> mGlobalWeights; // sum of peptide-length weighted counts from all alleles
    private final Map<Integer,Double> mGlobalPeptideLengthTotals;

    public PosWeightModel(final CalcConstants calcConstants, final HlaSequences hlaSequences)
    {
        mConstants = calcConstants;
        mAminoAcidFrequency = new AminoAcidFrequency();
        mBlosumMapping = new BlosumMapping();
        mHlaSequences = hlaSequences;
        mNoiseModel = new NoiseModel(mAminoAcidFrequency, calcConstants.NoiseProbability, calcConstants.NoiseWeight);
        mGlobalWeights = Maps.newHashMap();
        mGlobalPeptideLengthTotals = Maps.newHashMap();
    }

    public boolean noiseEnabled() { return mNoiseModel.enabled(); }
    public final NoiseModel noiseModel() { return mNoiseModel; }

    public static final int INVALID_POS = -1;

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

        if(mConstants.GlobalWeight > 0)
        {
            double[][] globalCounts = mGlobalWeights.get(bindCounts.PeptideLength);

            if(globalCounts == null)
            {
                globalCounts = new double[AMINO_ACID_COUNT][bindCounts.PeptideLength];
                mGlobalWeights.put(bindCounts.PeptideLength, globalCounts);
            }

            for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
            {
                for(int pos = 0; pos < bindCounts.PeptideLength; ++pos)
                {
                    globalCounts[aa][pos] += weightedCounts[aa][pos];
                }
            }
        }
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
                finalWeightedCounts[aa][pos] = weightedCounts[aa][pos];

                for(BindCountData otherBindCounts : allBindCounts)
                {
                    if(otherBindCounts.Allele.equals(bindCounts.Allele))
                        continue;

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

    public void setGlobalTotals()
    {
        for(Map.Entry<Integer,double[][]> entry : mGlobalWeights.entrySet())
        {
            int peptideLength = entry.getKey();
            double[][] counts = entry.getValue();

            double total = calcCountsTotal(counts);
            mGlobalPeptideLengthTotals.put(peptideLength, total);
        }
    }

    private static double calcCountsTotal(final double[][] counts)
    {
        double total = 0;

        for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
        {
            for(int pos = 0; pos < counts[0].length; ++pos)
            {
                total += counts[aa][pos];
            }
        }

        return total;
    }

    public BindScoreMatrix createMatrix(final BindCountData bindCounts)
    {
        double globalWeight = mConstants.GlobalWeight;

        if(mGlobalPeptideLengthTotals.isEmpty() && globalWeight > 0)
            setGlobalTotals();

        final double[][] finalWeightedCounts = bindCounts.getFinalWeightedCounts();

        BindScoreMatrix matrix = new BindScoreMatrix(bindCounts.Allele, bindCounts.PeptideLength);
        final double[][] data = matrix.getBindScores();

        double globalReductionFactor = 1;
        double[][] globalCounts = mGlobalWeights.get(bindCounts.PeptideLength);

        if(globalWeight > 0)
        {
            double total = calcCountsTotal(bindCounts.getFinalWeightedCounts());
            double globalTotal = mGlobalPeptideLengthTotals.get(bindCounts.PeptideLength);
            globalReductionFactor = total / globalTotal;
        }

        for(int pos = 0; pos < bindCounts.PeptideLength; ++pos)
        {
            double posTotalCount = 0;

            for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
            {
                posTotalCount += finalWeightedCounts[aa][pos];
            }

            for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
            {
                double aaFrequency = mAminoAcidFrequency.getAminoAcidFrequency(aa);

                // Peptide Weight = log(max(posWeight(aa,pos) * AaAdjust, 0.005 * posWeightTotal)/AaFreq,2)
                double adjustedCount = finalWeightedCounts[aa][pos];

                if(aa == 'C')
                    adjustedCount *= AMINO_ACID_C_FREQ_ADJUST;

                // handle very low observation counts
                adjustedCount = max(adjustedCount, posTotalCount * MIN_OBSERVED_AA_POS_FREQ);

                if(globalWeight > 0)
                {
                    adjustedCount = globalWeight * globalReductionFactor * globalCounts[aa][pos] + (1 - globalWeight) * adjustedCount;
                }

                double posWeight = log(2, adjustedCount / aaFrequency);
                data[aa][pos] = posWeight;
            }
        }

        return matrix;
    }

    public void writeGlobalCounts(final BufferedWriter writer, int maxPeptideLength)
    {
        try
        {
            for(Map.Entry<Integer,double[][]> entry : mGlobalWeights.entrySet())
            {
                int peptideLength = entry.getKey();
                double[][] counts = entry.getValue();

                for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
                {
                    char aminoAcid = AMINO_ACIDS.get(aa);

                    writer.write(String.format("%s,%s,%d,%c", "Global", "GLOBAL", peptideLength, aminoAcid));

                    for(int pos = 0; pos < maxPeptideLength; ++pos)
                    {
                        if(pos < peptideLength)
                        {
                            writer.write(String.format(",%.1f", counts[aa][pos]));
                        }
                        else
                        {
                            writer.write(",0.0");
                        }
                    }

                    writer.newLine();
                }
            }
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write global counts data: {}", e.toString());
        }
    }

}
