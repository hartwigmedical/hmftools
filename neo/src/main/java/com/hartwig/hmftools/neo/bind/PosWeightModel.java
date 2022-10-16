package com.hartwig.hmftools.neo.bind;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.VectorUtils.clear;
import static com.hartwig.hmftools.neo.bind.BindConstants.ALLELE_POS_MAPPING_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACID_COUNT;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACID_C_FREQ_ADJUST;
import static com.hartwig.hmftools.neo.bind.BindConstants.MIN_OBSERVED_AA_POS_FREQ;
import static com.hartwig.hmftools.neo.bind.BindConstants.REF_PEPTIDE_LEFT_FIXED_POS;
import static com.hartwig.hmftools.neo.bind.BindConstants.REF_PEPTIDE_LENGTH;

import static org.apache.commons.math3.util.FastMath.log;

import java.util.List;

import com.hartwig.hmftools.common.aminoacid.BlosumMapping;
import com.hartwig.hmftools.neo.utils.AminoAcidFrequency;

public class PosWeightModel
{
    private final CalcConstants mConstants;
    private final AminoAcidFrequency mAminoAcidFrequency;
    private final BlosumMapping mBlosumMapping;
    private final HlaSequences mHlaSequences;
    private final NoiseModel mNoiseModel;

    public PosWeightModel(final CalcConstants calcConstants, final HlaSequences hlaSequences)
    {
        mConstants = calcConstants;
        mAminoAcidFrequency = new AminoAcidFrequency();
        mBlosumMapping = new BlosumMapping();
        mHlaSequences = hlaSequences;
        mNoiseModel = new NoiseModel(mAminoAcidFrequency, calcConstants.NoiseProbability, calcConstants.NoiseWeight);
    }

    public boolean noiseEnabled() { return mNoiseModel.enabled(); }
    public final NoiseModel noiseModel() { return mNoiseModel; }

    public static final int INVALID_POS = -1;

    public static int peptidePosToRef(int peptideLength, int position)
    {
        return peptidePosToRef(REF_PEPTIDE_LENGTH, peptideLength, position);
    }

    public static int peptidePosToRef(int refLength, int peptideLength, int position)
    {
        if(position <= REF_PEPTIDE_LEFT_FIXED_POS)
            return position;

        // padded from the right
        return position + (refLength - peptideLength);
    }

    public static int refPeptidePosToActual(int peptideLength, int refPosition)
    {
        return refPeptidePosToActual(REF_PEPTIDE_LENGTH, peptideLength, refPosition);
    }

    public static int refPeptidePosToActual(int refLength, int peptideLength, int refPosition)
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

    public void buildWeightedCounts(final BindCountData bindCounts, final List<BindCountData> pepLenBindCounts)
    {
        // translate the counts from various peptide lengths into this set of counts, normalising to the reference peptide length

        // PWF = PeptideLengthFactor parameter
        // PWM = PeptideLengthMaxWeight parameter

        // RC = raw bind count per amino acid
        // AOC = adjusted other-length AA count = otherCount * PWF
        // AOTC = adjusted other-length total count = sum(all AAs: otherCount) * PWF

        // result: PLWC = peptide-length weighted count
        // PLWC = own-count + sum(AOC each PL) P13 * PLMW / max(PLMW, sum(AOTC each PL))

        final double[][] counts = mNoiseModel.enabled() ? bindCounts.getNoiseCounts() : bindCounts.getBindCounts();
        final double[][] weightedCounts = bindCounts.getWeightedCounts();

        double[] adjustedAACounts = new double[AMINO_ACID_COUNT];

        final double weightFactor = mConstants.PeptideLengthWeightFactor;
        final double weightMax = mConstants.PeptideLengthWeightMax;

        for(int pos = 0; pos < bindCounts.PeptideLength; ++pos)
        {
            // first get the totals across the other peptide-length counts
            double adjOtherTotal = 0;

            clear(adjustedAACounts);

            for(BindCountData otherBindCounts : pepLenBindCounts)
            {
                if(otherBindCounts.PeptideLength == bindCounts.PeptideLength)
                    continue;

                final double[][] otherCounts = mNoiseModel.enabled() ? otherBindCounts.getNoiseCounts() : otherBindCounts.getBindCounts();

                int refPos = peptidePosToRef(bindCounts.PeptideLength, pos);
                int otherPos = refPeptidePosToActual(otherBindCounts.PeptideLength, refPos);

                if(otherPos == INVALID_POS) // the other counts do not have this position, so ignore it for them
                    continue;

                for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
                {
                    double otherCount = otherCounts[aa][otherPos];
                    double adjOtherCount = otherCount * weightFactor;
                    adjustedAACounts[aa] += adjOtherCount;
                    adjOtherTotal += adjOtherCount;
                }
            }

            for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
            {
                weightedCounts[aa][pos] = counts[aa][pos] + weightFactor * adjustedAACounts[aa] * weightMax / max(weightMax, adjOtherTotal);
            }
        }
    }

    public void buildPositionAdjustedTotals(final List<BindCountData> pepLenBindCounts)
    {
        final double weightFactor = mConstants.PeptideLengthWeightFactor;
        final double weightMax = mConstants.PeptideLengthWeightMax;

        for(BindCountData bindCounts : pepLenBindCounts)
        {
            double[] pepLenPosTotals = new double[bindCounts.PeptideLength];
            double[] otherPepLenPosTotals = new double[bindCounts.PeptideLength];

            for(BindCountData otherBindCounts : pepLenBindCounts)
            {
                final double[][] otherCounts =
                        mNoiseModel.enabled() ? otherBindCounts.getNoiseCounts() : otherBindCounts.getBindCounts();

                if(otherBindCounts.PeptideLength == bindCounts.PeptideLength)
                {
                    for(int pos = 0; pos < bindCounts.PeptideLength; ++pos)
                    {
                        for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
                        {
                            pepLenPosTotals[pos] += otherCounts[aa][pos];
                        }
                    }
                }
                else
                {
                    for(int pos = 0; pos < bindCounts.PeptideLength; ++pos)
                    {
                        int refPos = peptidePosToRef(bindCounts.PeptideLength, pos);
                        int otherPos = refPeptidePosToActual(otherBindCounts.PeptideLength, refPos);

                        if(otherPos == INVALID_POS) // the other counts do not have this position, so ignore it for them
                            continue;

                        for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
                        {
                            double adjOtherCount = otherCounts[aa][otherPos] * weightFactor;
                            otherPepLenPosTotals[pos] += adjOtherCount;
                        }
                    }
                }
            }

            // adjusted pos-total = per peptide-length position total
            // adjPT = position total for peptide-length + min(PWM, PWF * other peptide-length position totals)
            final double[] adjustedPosTotals = bindCounts.getAdjustedPosTotals();

            for(int pos = 0; pos < bindCounts.PeptideLength; ++pos)
            {
                adjustedPosTotals[pos] = pepLenPosTotals[pos] + min(weightMax, weightFactor * otherPepLenPosTotals[pos]);
            }
        }
    }

    public void buildFinalWeightedCounts(final BindCountData bindCounts, final List<BindCountData> allBindCounts)
    {
        // calculate Blosum similarity at this position vs all the other alleles from matching peptide lengths

        // WCount(A,L,P,AA) = LWCount(A,L,P,AA) + MWF* SUM(a<>A)  [ LWCount(a,L,P,AA) * (2^(LogSim(m,M)) /  MAX(i=all motifs)[2^(LogSim(i,M))]]
        // * maxMW / max(MWF * SUM(a<>A)  [ LWCount(a,L,P) * (2^(LogSim(m,M)) /  MAX(i=all motifs)[2^(LogSim(i,M))]],maxMW)

        // PWF = PeptideLengthWeightFactor parameter
        // PWM = PeptideLengthWeightMax parameter
        // AMWF = AlleleMotifWeightFactor parameter
        // AMW = AlleleMotifWeightMax parameter

        // AMWC = own PLWC + sum(other adjPLWCs) * AMW / max(AMW, sum(adjPT * motifSimilarity - all alleles))

        // sumIFS(Q:Q,L:L,B1) + sumIFS(Q:Q,L:L,"<>"&B1) * B7/max(B7,sum(R:R))

        // where:
        // PLWC = peptide-length weighted count - as previously computed per AA and position
        // adjPLWC = AMWF * other-allele PLWC * motifSimilarity
        // adjPT as previously computed per allele, peptide length and position

        double alleleWeightFactor = mConstants.AlleleWeightFactor;
        double alleleWeightMax = mConstants.AlleleWeightMax;

        final double[][] weightedCounts = bindCounts.getWeightedCounts();
        final double[][] finalWeightedCounts = bindCounts.getFinalWeightedCounts();

        for(int pos = 0; pos < bindCounts.PeptideLength; ++pos)
        {
            // pos needs to be mapped to the 9-mer allele position peptide-position mapping
            int refPos = peptidePosToRef(bindCounts.PeptideLength, pos);
            int mappingPos = refPeptidePosToActual(ALLELE_POS_MAPPING_PEPTIDE_LENGTH, refPos);

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

            double selfScore = mBlosumMapping.calcSequenceScore(positionMotif);

            // cache motif similarities per position once since used repeatedly below
            final double[] motifSimilarities = new double[allBindCounts.size()];

            // calculate denominator
            double adjPosSimilarityTotal = 0;

            for(int alleleIndex = 0; alleleIndex < allBindCounts.size(); ++alleleIndex)
            {
                BindCountData otherBindCounts = allBindCounts.get(alleleIndex);
                String otherPositionMotif = mHlaSequences.getSequence(otherBindCounts.Allele, mappingPos);

                if(otherPositionMotif == null)
                    continue;

                double motifSimilarity = 1;

                if(!otherPositionMotif.equals(positionMotif))
                {
                    double crossAlleleScore = mBlosumMapping.calcSequenceScore(positionMotif, otherPositionMotif);
                    motifSimilarity = crossAlleleScore / selfScore;
                }

                motifSimilarities[alleleIndex] = motifSimilarity;

                if(otherBindCounts.Allele.equals(bindCounts.Allele))
                    adjPosSimilarityTotal += otherBindCounts.getAdjustedPosTotals()[pos];
                else
                    adjPosSimilarityTotal += otherBindCounts.getAdjustedPosTotals()[pos] * motifSimilarity * alleleWeightFactor;
            }

            double adjPosSimilarityTotalFactor = alleleWeightMax / max(alleleWeightMax, adjPosSimilarityTotal);

            for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
            {
                for(int alleleIndex = 0; alleleIndex < allBindCounts.size(); ++alleleIndex)
                {
                    BindCountData otherBindCounts = allBindCounts.get(alleleIndex);

                    if(otherBindCounts.Allele.equals(bindCounts.Allele))
                    {
                        // taken as is, no weighting
                        finalWeightedCounts[aa][pos] += weightedCounts[aa][pos];
                        continue;
                    }

                    double otherCount = otherBindCounts.getWeightedCounts()[aa][pos];

                    if(otherCount == 0)
                        continue;

                    double motifSimilarity = motifSimilarities[alleleIndex];

                    double otherWeightedCount = otherCount * motifSimilarity * alleleWeightFactor * adjPosSimilarityTotalFactor;
                    finalWeightedCounts[aa][pos] += otherWeightedCount;
                }
            }
        }
    }

    public BindScoreMatrix createMatrix(final BindCountData bindCounts)
    {
        final double[][] finalWeightedCounts = bindCounts.getFinalWeightedCounts();

        BindScoreMatrix matrix = new BindScoreMatrix(bindCounts.Allele, bindCounts.PeptideLength);
        final double[][] data = matrix.getBindScores();

        for(int pos = 0; pos < bindCounts.PeptideLength; ++pos)
        {
            double posTotalCount = 0;

            for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
            {
                posTotalCount += finalWeightedCounts[aa][pos];
            }

            posTotalCount = max(posTotalCount, 0.001); // for training data without observed counts in a peptide length

            for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
            {
                double aaFrequency = mAminoAcidFrequency.getAminoAcidFrequency(aa);

                // Peptide Weight = log(max(posWeight(aa,pos), 0.005 * posWeightTotal)/AaFreq,2)
                // Peptide Weight new = log(max(posWeight(aa,pos)/posWeightTotal, 0.005)/AaFreq,2) - now a % weight
                double adjustedCount = finalWeightedCounts[aa][pos];

                if(aa == 'C')
                    adjustedCount *= AMINO_ACID_C_FREQ_ADJUST;

                // handle very low observation counts
                adjustedCount = max(adjustedCount / posTotalCount, MIN_OBSERVED_AA_POS_FREQ);

                double posWeight = log(2, adjustedCount / aaFrequency);
                data[aa][pos] = posWeight;
            }
        }

        return matrix;
    }

}
