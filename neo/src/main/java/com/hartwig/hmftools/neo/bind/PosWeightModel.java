package com.hartwig.hmftools.neo.bind;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.pow;

import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.ALLELE_POS_MAPPING_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACIDS;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACID_COUNT;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACID_C_FREQ_ADJUST;
import static com.hartwig.hmftools.neo.bind.BindConstants.MIN_OBSERVED_AA_POS_FREQ;
import static com.hartwig.hmftools.neo.bind.BindConstants.REF_PEPTIDE_LEFT_FIXED_POS;
import static com.hartwig.hmftools.neo.bind.BindConstants.REF_PEPTIDE_LENGTH;

import static org.apache.commons.math3.util.FastMath.log;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.neo.utils.AminoAcidFrequency;

public class PosWeightModel
{
    private final CalcConstants mConstants;
    private final AminoAcidFrequency mAminoAcidFrequency;
    private final BlosumMapping mBlosumMapping;
    private final HlaSequences mHlaSequences;

    public PosWeightModel(final CalcConstants calcConstants, final HlaSequences hlaSequences)
    {
        mConstants = calcConstants;
        mAminoAcidFrequency = new AminoAcidFrequency();
        mBlosumMapping = new BlosumMapping();
        mHlaSequences = hlaSequences;
    }

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

    public void buildWeightedCounts(final BindCountData bindCounts, final List<BindCountData> peptideLengthCounts)
    {
        // translate the counts from various peptide lengths into this set of counts, normalising to the reference peptide length

        final double[][] weightedCounts = bindCounts.getWeightedCounts();

        for(BindCountData otherBindCounts : peptideLengthCounts)
        {
            final double[][] otherCounts = otherBindCounts.getBindCounts();
            int otherPeptideLength = otherBindCounts.PeptideLength;

            // LWCount(A,L,P,AA) = Count(A,L,P,AA) + SUM(l<>L) [ Count(A,l,P,AA) * 1/abs(L-l)] * [1 / ( 1+ Obs(A,L)/LHW)^E)]
            double weight = 1;

            if(otherPeptideLength != bindCounts.PeptideLength)
            {
                weight = 1.0 / abs(bindCounts.PeptideLength - otherPeptideLength)
                        * 1 / pow(1 + bindCounts.totalBindCount() / mConstants.PeptideLengthWeight, mConstants.WeightExponent);
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
    }

    public void buildFinalWeightedCounts(
            final BindCountData bindCounts, final List<BindCountData> allBindCounts, final Map<String,Integer> alleleTotalCounts)
    {
        // calculate blosum similarity at this position vs all the other alleles from matching peptide lengths

        // WCount(A,L,P,AA) = LWCount(A,L,P,AA) + SUM(a<>A)  [ LWCount(a,L,P,AA)
        // * (2^(LogSim(m,M)) /  MAX(i=all motifs)[2^(LogSim(i,M))]] * [1 / ( 1+ Obs(A,L)/MHW)^E)]

        final double[][] weightedCounts = bindCounts.getWeightedCounts();
        final double[][] finalWeightedCounts = bindCounts.getFinalWeightedCounts();

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

                    int alleleCount = alleleTotalCounts.get(otherBindCounts.Allele);
                    double observationsWeight = 1 / pow(1 + alleleCount / mConstants.AlleleWeight, mConstants.WeightExponent);

                    double otherWeightedCount = otherCount * observationsWeight * motifSimilarity;
                    finalWeightedCounts[aa][pos] += otherWeightedCount;
                }
            }
        }
    }

    public BindScoreMatrix createMatrix(final BindCountData bindCounts, final Map<Integer,Integer> peptideLengthFrequency)
    {
        NE_LOGGER.debug("creating allele({}) peptideLength({}) matrix data", bindCounts.Allele, bindCounts.PeptideLength);

        final double[][] finalWeightedCounts = bindCounts.getFinalWeightedCounts();

        BindScoreMatrix matrix = new BindScoreMatrix(bindCounts.Allele, bindCounts.PeptideLength);
        final double[][] data = matrix.getBindScores();

        double alleleBinds = peptideLengthFrequency.values().stream().mapToInt(x -> x.intValue()).sum();
        double alleleLengthPerc = bindCounts.totalBindCount() / alleleBinds;

        for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
        {
            char aminoAcid = AMINO_ACIDS.get(aa);
            double aaFrequency = mAminoAcidFrequency.getAminoAcidFrequency(aminoAcid);

            for(int pos = 0; pos < bindCounts.PeptideLength; ++pos)
            {
                // Peptide Weight = 1/L * Sum[Log2(P(x,i)/Q(x) * Obs(L,A) / SUM[Obs(l,A)])]
                double adjustedCount = finalWeightedCounts[aa][pos];

                if(aa == 'C')
                    adjustedCount *= AMINO_ACID_C_FREQ_ADJUST;

                // handle very low observation counts
                adjustedCount = max(adjustedCount, bindCounts.totalBindCount() * MIN_OBSERVED_AA_POS_FREQ);

                double weightedCount = log(2, adjustedCount / aaFrequency);
                double posWeight = weightedCount; //  * (1.0 / PeptideLength); // * (mTotalBinds / alleleBinds);
                data[aa][pos] = posWeight;
            }
        }

        return matrix;
    }

}
