package com.hartwig.hmftools.neo.bind;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.pow;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.ALLELE_POS_MAPPING_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACIDS;
import static com.hartwig.hmftools.neo.bind.BindConstants.INVALID_AMINO_ACID;
import static com.hartwig.hmftools.neo.bind.BindConstants.MIN_OBSERVED_AA_POS_FREQ;
import static com.hartwig.hmftools.neo.bind.BindConstants.REF_PEPTIDE_LEFT_FIXED_POS;
import static com.hartwig.hmftools.neo.bind.BindConstants.REF_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.BindConstants.aminoAcidIndex;

import static org.apache.commons.math3.util.FastMath.log;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.neo.utils.AminoAcidFrequency;

public class BindCountData
{
    public final String Allele;
    public final int PeptideLength;

    private final int mAminoAcidCount;

    private final int[][] mObservations; // data points with a given amino acid and peptide
    private final double[][] mBindScoreTotals; // log-score total
    private final double[][] mBindCounts; // observations with affinity below configured binding threshold

    private final double[][] mWeightedCounts; // weighted across the peptide lengths for this allele, padded to reference length
    private final double[][] mFinalWeightedCounts; // weighted across all other alleles, padded to reference length

    private final Map<String,double[]> mComboData;
    private int mTotal;
    private int mTotalBinds; // vs high affinity threshold
    private double mCalcTotalBinds; // vs high affinity threshold

    public BindCountData(final String allele, final int peptideLength)
    {
        Allele = allele;
        PeptideLength = peptideLength;
        mAminoAcidCount = AMINO_ACIDS.size();

        mObservations = new int[mAminoAcidCount][PeptideLength];
        mBindScoreTotals = new double[mAminoAcidCount][PeptideLength];
        mBindCounts = new double[mAminoAcidCount][PeptideLength];
        mWeightedCounts = new double[mAminoAcidCount][PeptideLength];
        mFinalWeightedCounts = new double[mAminoAcidCount][PeptideLength];

        mComboData = Maps.newHashMap();
        mTotal = 0;
        mTotalBinds = 0;
        mCalcTotalBinds = 0;
    }

    public final double[][] getBindCounts() { return mBindCounts; }
    public final int[][] getObservations() { return mObservations; }
    public Map<String,double[]> getComboData() { return mComboData; }
    public final double[][] getWeightedCounts() { return mWeightedCounts; }
    public final double[][] getFinalWeightedCounts() { return mFinalWeightedCounts; }
    public int totalBindCount() { return mTotalBinds; }

    public void logStats()
    {
        NE_LOGGER.info("allele({}) peptideLen({}) total({}) binds({}) pairs({})",
                Allele, PeptideLength, mTotal, mTotalBinds, mComboData.size());
    }

    public void processBindData(final BindData bindData, boolean calcPairs, final CalcConstants calcConstants)
    {
        if(bindData.Peptide.length() != PeptideLength || !bindData.Allele.equals(Allele))
            return;

        double levelScore = calcConstants.deriveLevelScore(bindData.Affinity);
        double bindPerc = calcConstants.deriveAffinityPercent(bindData.Affinity);

        ++mTotal;
        mCalcTotalBinds += bindPerc;

        boolean actualBind = bindData.Affinity < calcConstants.BindingAffinityHigh;
        boolean predictedBind = bindData.predictedAffinity() < calcConstants.BindingAffinityHigh && actualBind;

        if(actualBind)
            ++mTotalBinds;

        for(int pos = 0; pos < bindData.Peptide.length(); ++pos)
        {
            char aminoAcid = bindData.Peptide.charAt(pos);
            int aaIndex = aminoAcidIndex(aminoAcid);

            if(aaIndex == INVALID_AMINO_ACID)
                continue;

            ++mObservations[aaIndex][pos];
            mBindScoreTotals[aaIndex][pos] += levelScore;
            mBindCounts[aaIndex][pos] += bindPerc;

            if(calcPairs && bindData.isTraining())
            {
                if(pos < bindData.Peptide.length() - 1)
                {
                    for(int pos2 = pos + 1; pos2 < bindData.Peptide.length(); ++pos2)
                    {
                        char aminoAcid2 = bindData.Peptide.charAt(pos2);
                        if(aminoAcidIndex(aminoAcid2) == INVALID_AMINO_ACID)
                            continue;

                        ComboCorrelations.updatePairData(
                                mComboData, aminoAcid, pos, aminoAcid2, pos2, levelScore, actualBind, predictedBind);
                    }
                }
            }
        }
    }

    public void processBindingPeptide(final String peptide)
    {
        // simpler method assumes the correct allele and that peptide binds
        ++mTotal;
        ++mTotalBinds;

        for(int pos = 0; pos < peptide.length(); ++pos)
        {
            char aminoAcid = peptide.charAt(pos);
            int aaIndex = aminoAcidIndex(aminoAcid);

            if(aaIndex == INVALID_AMINO_ACID)
                continue;

            ++mObservations[aaIndex][pos];
            mBindScoreTotals[aaIndex][pos] += 1;
            mBindCounts[aaIndex][pos] += 1;
        }
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

    public void buildWeightedCounts(final List<BindCountData> peptideLengthCounts, final CalcConstants calcConst)
    {
        // translate the counts from various peptide lengths into this set of counts, normalising to the reference peptide length
        for(BindCountData otherBindCounts : peptideLengthCounts)
        {
            final double[][] otherCounts = otherBindCounts.getBindCounts();
            int otherPeptideLength = otherBindCounts.PeptideLength;

            // LWCount(A,L,P,AA) = Count(A,L,P,AA) + SUM(l<>L) [ Count(A,l,P,AA) * 1/abs(L-l)] * [1 / ( 1+ Obs(A,L)/LHW)^E)]
            double weight = 1;

            if(otherPeptideLength != PeptideLength)
            {
                weight = 1.0 / abs(PeptideLength - otherPeptideLength)
                        * 1 / pow(1 + mTotalBinds / calcConst.PeptideLengthWeight, calcConst.WeightExponent);
            }

            for(int aa = 0; aa < mAminoAcidCount; ++aa)
            {
                for(int otherPos = 0; otherPos < otherPeptideLength; ++otherPos)
                {
                    int refPos = peptidePositionToRef(REF_PEPTIDE_LENGTH, otherPeptideLength, otherPos);
                    int pos = refPeptidePositionToActual(REF_PEPTIDE_LENGTH, PeptideLength, refPos);

                    if(pos == INVALID_POS)
                        continue;

                    mWeightedCounts[aa][pos] += weight * otherCounts[aa][otherPos];
                }
            }
        }
    }

    public void buildFinalWeightedCounts(
            final List<BindCountData> allBindCounts, final Map<String,Integer> alleleTotalCounts, final CalcConstants calcConst,
            final BlosumMapping blosumMapping, final HlaSequences hlaSequences)
    {
        // calculate blosum similarity at this position vs all the other alleles from matching peptide lengths

        // WCount(A,L,P,AA) = LWCount(A,L,P,AA) + SUM(a<>A)  [ LWCount(a,L,P,AA)
        // * (2^(LogSim(m,M)) /  MAX(i=all motifs)[2^(LogSim(i,M))]] * [1 / ( 1+ Obs(A,L)/MHW)^E)]

        for(int pos = 0; pos < PeptideLength; ++pos)
        {
            // pos needs to be mapped to the 9-mer allele position peptide-position mapping
            int refPos = peptidePositionToRef(REF_PEPTIDE_LENGTH, PeptideLength, pos);
            int mappingPos = refPeptidePositionToActual(REF_PEPTIDE_LENGTH, ALLELE_POS_MAPPING_PEPTIDE_LENGTH, refPos);

            if(mappingPos == INVALID_POS)
            {
                // cannot use other alleles
                for(int aa = 0; aa < mAminoAcidCount; ++aa)
                {
                    mFinalWeightedCounts[aa][pos] = mWeightedCounts[aa][pos];
                }

                continue;
            }

            String positionMotif = hlaSequences.getSequence(Allele, mappingPos);

            if(positionMotif == null)
                continue;

            double selfScore = blosumMapping.calcSequenceBlosumScore(positionMotif);

            for(int aa = 0; aa < mAminoAcidCount; ++aa)
            {
                mFinalWeightedCounts[aa][pos] = mWeightedCounts[aa][pos];

                for(BindCountData otherBindCounts : allBindCounts)
                {
                    if(otherBindCounts.Allele.equals(Allele))
                        continue;

                    double otherCount = otherBindCounts.getWeightedCounts()[aa][pos];

                    if(otherCount == 0)
                        continue;

                    String otherPositionMotif = hlaSequences.getSequence(otherBindCounts.Allele, mappingPos);

                    if(otherPositionMotif == null)
                        continue;

                    double motifSimilarity = 1;

                    if(!otherPositionMotif.equals(positionMotif))
                    {
                        double crossAlleleScore = blosumMapping.calcSequenceBlosumScore(positionMotif, otherPositionMotif);
                        motifSimilarity = crossAlleleScore / selfScore;
                    }

                    int alleleCount = alleleTotalCounts.get(otherBindCounts.Allele);
                    double observationsWeight = 1 / pow(1 + alleleCount / calcConst.AlleleWeight, calcConst.WeightExponent);

                    double otherWeightedCount = otherCount * observationsWeight * motifSimilarity;
                    mFinalWeightedCounts[aa][pos] += otherWeightedCount;
                }
            }
        }
    }

    public BindScoreMatrix createMatrix(final AminoAcidFrequency aminoAcidFrequency, final Map<Integer,Integer> peptideLengthFrequency)
    {
        NE_LOGGER.debug("creating allele({}) peptideLength({}) matrix data", Allele, PeptideLength);

        BindScoreMatrix matrix = new BindScoreMatrix(Allele, PeptideLength);
        final double[][] data = matrix.getBindScores();

        double alleleBinds = peptideLengthFrequency.values().stream().mapToInt(x -> x.intValue()).sum();
        double alleleLengthPerc = mTotalBinds / alleleBinds;

        for(int aa = 0; aa < mAminoAcidCount; ++aa)
        {
            char aminoAcid = AMINO_ACIDS.get(aa);
            double aaFrequency = aminoAcidFrequency.getAminoAcidFrequency(aminoAcid);

            for(int pos = 0; pos < PeptideLength; ++pos)
            {
                // Peptide Weight = 1/L * Sum[Log2(P(x,i)/Q(x) * Obs(L,A) / SUM[Obs(l,A)])]

                double adjustedCount = max(mFinalWeightedCounts[aa][pos], mTotalBinds * MIN_OBSERVED_AA_POS_FREQ);
                double weightedCount = log(2, adjustedCount / aaFrequency);
                double posWeight = weightedCount; //  * (1.0 / PeptideLength); // * (mTotalBinds / alleleBinds);
                data[aa][pos] = posWeight;

                /* simple original scoring
                double totalBinds = mBindCounts[aa][pos];
                double freqPerc = totalBinds / mCalcTotalBinds;
                freqPerc = max(freqPerc, MIN_OBSERVED_AA_POS_FREQ);

                // Score = log(max(FreqPerc,0.001) / AminoAcidFreq,2)
                double score = log(2, freqPerc / aaFrequency);

                data[aa][pos] = score;
                */
            }
        }

        return matrix;
    }

    public static BufferedWriter initMatrixWriter(final String filename, int peptideLength, boolean writeCounts)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Allele,PeptideLength,AminoAcid");

            if(writeCounts)
                writer.write(",DataType");

            for(int i = 0; i < peptideLength; ++i)
            {
                writer.write(String.format(",P%d", i));
            }

            writer.newLine();

            return writer;
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to initialise matrix data file({}): {}", filename, e.toString());
            return null;
        }
    }

    public void writeMatrixData(final BufferedWriter writer, final BindScoreMatrix matrix, int maxPeptideLength, boolean writeCounts)
    {
        NE_LOGGER.debug("writing allele({}) matrix data", Allele);

        final double[][] data = matrix.getBindScores();

        try
        {
            for(int aa = 0; aa < mAminoAcidCount; ++aa)
            {
                char aminoAcid = AMINO_ACIDS.get(aa);

                writer.write(String.format("%s,%d,%c", Allele, PeptideLength, aminoAcid));

                if(writeCounts)
                    writer.write(",PosWeights");

                for(int pos = 0; pos < maxPeptideLength; ++pos)
                {
                    if(pos < PeptideLength)
                    {
                        writer.write(String.format(",%.6f", data[aa][pos]));
                    }
                    else
                    {
                        writer.write(",0.0");
                    }
                }

                writer.newLine();
            }

            if(writeCounts)
            {
                for(int i = 0; i <= 2; ++i)
                {
                    final double[][] counts = (i == 0) ? mBindCounts : (i == 1) ? mWeightedCounts : mFinalWeightedCounts;
                    String dataType = (i == 0) ? "BindCounts" : (i == 1) ? "PeptideLengthWeighted" : "AlleleMotifWeighted";

                    for(int aa = 0; aa < mAminoAcidCount; ++aa)
                    {
                        char aminoAcid = AMINO_ACIDS.get(aa);

                        writer.write(String.format("%s,%d,%c,%s", Allele, PeptideLength, aminoAcid, dataType));

                        for(int pos = 0; pos < maxPeptideLength; ++pos)
                        {
                            if(pos < PeptideLength)
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
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write matrix data: {}", e.toString());
        }
    }

    public static BufferedWriter initFrequencyWriter(final String filename)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Allele,PeptideLength,AminoAcid,PeptidePos,Count,Score,TotalScore,ActualBinds");
            writer.newLine();

            return writer;
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to initialise frequency data file({}): {}", filename, e.toString());
            return null;
        }
    }

    public void writeFrequencyData(final BufferedWriter writer)
    {
        NE_LOGGER.debug("writing allele({}) frequency", Allele);

        try
        {
            for(int aa = 0; aa < mAminoAcidCount; ++aa)
            {
                char aminoAcid = AMINO_ACIDS.get(aa);

                for(int p = 0; p < PeptideLength; ++p)
                {
                    int freq = mObservations[aa][p];
                    double totalBinds = mBindCounts[aa][p];
                    double totalScore = mBindScoreTotals[aa][p];
                    double avgScore = freq > 0 ? totalScore / freq : 0;

                    writer.write(String.format("%s,%d,%s,%d,%d,%.4f,%.4f,%.2f",
                            Allele, PeptideLength, aminoAcid, p, freq, avgScore, totalScore, totalBinds));

                    writer.newLine();
                }
            }
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write frequency data: {}", e.toString());
        }
    }

    public String toString()
    {
        return String.format("%s L%d: totalBinds(%d)", Allele, PeptideLength, mTotalBinds);
    }

}
