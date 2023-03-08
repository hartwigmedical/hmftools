package com.hartwig.hmftools.neo.scorer;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.codon.AminoAcidRna.AA_SELENOCYSTEINE;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.MIN_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.BindConstants.REF_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.FlankCounts.FLANK_AA_COUNT;
import static com.hartwig.hmftools.neo.scorer.NeoRnaData.NO_TPM_VALUE;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.neo.NeoEpitopeFile;
import com.hartwig.hmftools.common.neo.NeoEpitopeType;
import com.hartwig.hmftools.neo.PeptideData;
import com.hartwig.hmftools.neo.epitope.EpitopeUtils;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class TpmCalculator
{
    private final int mFlankCount;
    private final int[] mPeptideLengthRange;
    private final Map<Integer,PoissonRangeValues> mPoissonRangeValues;

    public static final double LOW_PROBABILITY = 0.05;
    public static final double HIGH_PROBABILITY = 0.95;

    private static final double EFFECTIVE_TPM_ACTUAL_PERC = 0.8;

    public static final int[] DEFAULT_PEPTIDE_LENGTH_RANGE = new int[] { MIN_PEPTIDE_LENGTH, REF_PEPTIDE_LENGTH };

    public TpmCalculator(final int flankCount, final int[] peptideLengthRange)
    {
        mFlankCount = flankCount;
        mPeptideLengthRange = peptideLengthRange;
        mPoissonRangeValues = Maps.newHashMap();
        populatePoissonCache();
    }

    public void compute(final String sampleId, final List<NeoEpitopeData> neoDataList)
    {
        double tpmNormalisationFactor = calculateTpmNormalisation(neoDataList);

        // derive the set of peptides per allele from the novel amino acids
        Map<String,List<NeoEpitopeData>> variantNeoMap = Maps.newHashMap();

        for(NeoEpitopeData neoData : neoDataList)
        {
            List<NeoEpitopeData> variantNeos = variantNeoMap.get(neoData.VariantInfo);

            if(variantNeos == null)
                variantNeoMap.put(neoData.VariantInfo, Lists.newArrayList(neoData));
            else
                variantNeos.add(neoData);
        }

        int peptideAlleleCount = 0;

        // evaluate peptides per variants to find repeated peptides from different transcripts
        for(List<NeoEpitopeData> variantNeos : variantNeoMap.values())
        {
            List<PeptideScoreData> peptideScoreDataList = Lists.newArrayList();

            double totalTpmDown = variantNeos.stream().mapToDouble(x -> x.RnaData.getTPM(FS_DOWN)).sum();

            for(NeoEpitopeData neoData : variantNeos)
            {
                double tpmUp = neoData.RnaData.getTPM(FS_UP);
                double tpmDown = neoData.RnaData.getTPM(FS_DOWN);

                List<PeptideData> peptides = EpitopeUtils.generatePeptides(
                        neoData.UpAminoAcids, neoData.NovelAminoAcids, neoData.DownAminoAcids, mPeptideLengthRange, mFlankCount);

                for(PeptideData peptide : peptides)
                {
                    if(peptide.Peptide.contains(AA_SELENOCYSTEINE))
                        continue;

                    PeptideScoreData peptideScoreData = peptideScoreDataList.stream()
                            .filter(x -> x.Peptide.equals(peptide.Peptide)).findFirst().orElse(null);

                    if(peptideScoreData != null)
                    {
                        peptideScoreData.addNeoepitopeTpm(neoData.Id, tpmDown);
                    }
                    else
                    {
                        peptideScoreDataList.add(new PeptideScoreData(peptide, neoData.Id, tpmUp, tpmDown));
                    }
                }
            }

            peptideAlleleCount += peptideScoreDataList.size();

            // now that the search for common peptides has been done, calculate expected and effective TPM values
            for(NeoEpitopeData neoData : variantNeos)
            {
                List<PeptideScoreData> neoPeptides = peptideScoreDataList.stream().filter(x -> x.hasNeo(neoData.Id)).collect(Collectors.toList());

                neoData.addPeptides(neoPeptides);

                double rawEffectiveTpm = tpmNormalisationFactor != NO_TPM_VALUE ?
                        neoData.RnaData.fragmentSupport() * tpmNormalisationFactor : NO_TPM_VALUE;

                PoissonRangeValues rangeValues = getOrCalcRangeValues(neoData.RnaData.fragmentSupport());

                for(PeptideScoreData peptideData : neoPeptides)
                {
                    // expected TPM is the sum of TPM Down from neoepitopes with this peptide / total TPM down for all neoepitopes of this variant
                    double expectedTpm = totalTpmDown > 0 ? peptideData.tpmUp() * peptideData.tpmDownTotal() / totalTpmDown : 0;

                    double effectiveTpm;

                    if(tpmNormalisationFactor == NO_TPM_VALUE || !neoData.RnaData.hasExpression())
                    {
                        effectiveTpm = expectedTpm;
                    }
                    else
                    {
                        // effectiveTPM = 0.8 * rawEffectiveTPM + 0.2 * max[poisson5%(rawEffectiveTPM), min[Poisson95%(rawEffectiveTPM), expectedTPM]]
                        double rangeLow = rangeValues.LowValue * tpmNormalisationFactor;
                        double rangeHigh = rangeValues.HighValue * tpmNormalisationFactor;

                        effectiveTpm = EFFECTIVE_TPM_ACTUAL_PERC * rawEffectiveTpm
                                + (1 - EFFECTIVE_TPM_ACTUAL_PERC) * max(rangeHigh, min(rangeLow, expectedTpm));
                    }

                    peptideData.setCalculatedTpms(rawEffectiveTpm, effectiveTpm, expectedTpm);
                }
            }
        }

        NE_LOGGER.debug("sample({}) neoepitopes({}) derived {} allele-peptides",
                sampleId, neoDataList.size(), peptideAlleleCount);
    }

    public static double calculateTpmNormalisation(final List<NeoEpitopeData> neoDataList)
    {
        double tpmUpTotal = 0;
        double baseDepthTotal = 0;
        int dataCount = 0;

        for(NeoEpitopeData neoData : neoDataList)
        {
            if(neoData.VariantType != NeoEpitopeType.MISSENSE)
                continue;

            if(neoData.RnaData.hasCoverage() && neoData.RnaData.hasExpression())
            {
                tpmUpTotal += neoData.RnaData.transExpression()[FS_UP];
                baseDepthTotal += neoData.RnaData.averageBaseDepth();
                ++dataCount;
            }
        }

        double tpmNormalisationFactor = dataCount > 0 ? tpmUpTotal / baseDepthTotal : NO_TPM_VALUE;
        return tpmNormalisationFactor;
    }

    private PoissonRangeValues getOrCalcRangeValues(double fragmentCount)
    {
        int fragCountClean = (int)round(fragmentCount);

        PoissonRangeValues rangeValues = mPoissonRangeValues.get(fragCountClean);

        if(rangeValues == null)
        {
            double lowValue = calcPoissonMean(fragCountClean, LOW_PROBABILITY);
            double highValue = calcPoissonMean(fragCountClean, HIGH_PROBABILITY);
            rangeValues = new PoissonRangeValues(lowValue, highValue);
            mPoissonRangeValues.put(fragCountClean, rangeValues);
        }

        return rangeValues;
    }

    private class PoissonRangeValues
    {
        public final double LowValue;
        public final double HighValue;

        public PoissonRangeValues(final double lowValue, final double highValue)
        {
            LowValue = lowValue;
            HighValue = highValue;
        }

        public String toString() { return format("low=%.2f high=%.2f", LowValue, HighValue); }
    }

    private void populatePoissonCache()
    {
        List<String> lines = new BufferedReader(new InputStreamReader(
                    TpmCalculator.class.getResourceAsStream("/ref/expected_fragment_ranges.csv"))).lines().collect(Collectors.toList());

        lines.remove(0);

        for(String line : lines)
        {
            String[] values = line.split(NeoEpitopeFile.DELIMITER, 3);

            mPoissonRangeValues.put(
                    Integer.parseInt( values[0]), new PoissonRangeValues(Double.parseDouble(values[1]), Double.parseDouble(values[2])));
        }
    }

    public static double calcPoissonMean(int fragmentCount, double requiredProb)
    {
        if(fragmentCount < 0)
            return 0;

        if(fragmentCount == 0)
            return requiredProb == HIGH_PROBABILITY ? 0 : 3;

        // find the mean for a Poisson distribution where the observed fragment count is at the required probability level
        int maxIterations = 20;
        int iterations = 0;

        double currentValue = 0;
        double testValueUpper = 0;
        double testValueLower = 0;

        if(requiredProb > 0.5)
        {
            currentValue = fragmentCount * 0.5;
            testValueUpper = fragmentCount;
            testValueLower = currentValue * 0.5;
        }
        else
        {
            currentValue = fragmentCount * 2;
            testValueUpper = currentValue * 3;
            testValueLower = fragmentCount;
        }

        PoissonDistribution poisson = new PoissonDistribution(currentValue);

        double currentProb = poisson.cumulativeProbability(fragmentCount);
        double probDiff = 0;

        while(iterations < maxIterations)
        {
            probDiff = abs(requiredProb - currentProb) / requiredProb;

            if(probDiff < 0.005)
                return currentValue;

            // if prob is too low, need to lower the test value
            if(currentProb < requiredProb)
            {
                testValueUpper = currentValue;
                currentValue = (currentValue + testValueLower) * 0.5;
            }
            else
            {
                testValueLower = currentValue;
                currentValue = (currentValue + testValueUpper) * 0.5;
            }

            poisson = new PoissonDistribution(currentValue);
            currentProb = poisson.cumulativeProbability(fragmentCount);
            ++iterations;
        }

        if(iterations >= maxIterations)
        {
            NE_LOGGER.warn(format("max iterations reached: value(%d) test(%d) prob(%.4f diff=%.4f)",
                    fragmentCount, currentValue, currentProb, probDiff));
        }

        return currentValue;
    }
}
