package com.hartwig.hmftools.neo.scorer;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.Math.sqrt;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.codon.AminoAcidRna.AA_SELENOCYSTEINE;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.MIN_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.BindConstants.REF_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.FlankCounts.FLANK_AA_COUNT;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.neo.NeoEpitopeType;
import com.hartwig.hmftools.neo.PeptideData;
import com.hartwig.hmftools.neo.epitope.EpitopeUtils;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class TpmCalculator
{
    private final Map<Double,PoissonRangeValues> mPoissonRangeValues;

    private static final int[] PEPTIDE_LENGTH_RANGE = new int[] { MIN_PEPTIDE_LENGTH, REF_PEPTIDE_LENGTH };

    private static final double LOW_PROBABILITY = 0.05;
    private static final double HIGH_PROBABILITY = 0.95;

    private static final double EFFECTIVE_TPM_ACTUAL_PERC = 0.8;

    public TpmCalculator()
    {
        mPoissonRangeValues = Maps.newHashMap();
    }

    public void compute(final String sampleId, final List<NeoEpitopeData> neoDataList, final Set<String> uniqueAlleles)
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

            for(NeoEpitopeData neoData : variantNeos)
            {
                double tpmUp = neoData.RnaData.getTPM(FS_UP);
                double tpmDown = neoData.RnaData.getTPM(FS_DOWN);

                List<PeptideData> peptides = EpitopeUtils.generatePeptides(
                        neoData.UpAminoAcids, neoData.NovelAminoAcids, neoData.DownAminoAcids, PEPTIDE_LENGTH_RANGE, FLANK_AA_COUNT);

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

                double rawEffectiveTpm = neoData.RnaData.fragmentSupport() * tpmNormalisationFactor;

                PoissonRangeValues rangeValues = getOrCalcRangeValues(rawEffectiveTpm);

                for(PeptideScoreData peptideData : neoPeptides)
                {
                    double expectedTpm = peptideData.tpmUpAllocation(neoData.Id);

                    // effectiveTPM = 0.8 * rawEffectiveTPM + 0.2 * max[poisson5%(rawEffectiveTPM), min[Poisson95%(rawEffectiveTPM), expectedTPM]]
                    double effectiveTpm = EFFECTIVE_TPM_ACTUAL_PERC * rawEffectiveTpm
                            + (1 - EFFECTIVE_TPM_ACTUAL_PERC) * max(rangeValues.LowValue, min(rangeValues.HighValue, expectedTpm));

                    peptideData.setEffectiveTpms(rawEffectiveTpm, effectiveTpm);
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

        double tpmNormalisationFactor = dataCount > 0 ? tpmUpTotal / baseDepthTotal : 1;
        return tpmNormalisationFactor;
    }

    private PoissonRangeValues getOrCalcRangeValues(double fragmentCount)
    {
        double fragCountClean = fragmentCount <= 5 ? round(fragmentCount/0.1) * 0.1 : round(fragmentCount);

        PoissonRangeValues rangeValues = mPoissonRangeValues.get(fragCountClean);

        if(rangeValues == null)
        {
            double lowValue = calcPoissonObservedGivenProb(fragCountClean, LOW_PROBABILITY);
            double highValue = calcPoissonObservedGivenProb(fragCountClean, HIGH_PROBABILITY);
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

        public String toString() { return format("low=%.1f high=%.1f", LowValue, HighValue); }
    }

    public static double calcPoissonObservedGivenProb(double expectedVal, double requiredProb)
    {
        if(expectedVal <= 0)
            return 0;

        if(expectedVal <= 1 && requiredProb == LOW_PROBABILITY)
            return 0;

        PoissonDistribution poisson = new PoissonDistribution(expectedVal);

        int maxIterations = 20;
        int iterations = 0;

        double refCount = 25;
        double refRangePerc = 0.44;
        double rangePerc = refRangePerc / sqrt(expectedVal/refCount);
        double range = rangePerc * expectedVal;

        int testValueUpper;
        int testValueLower;

        if(requiredProb > 0.5)
        {
            testValueUpper = (int) round(expectedVal + range * 1.5);
            testValueLower = (int) round(max(expectedVal + range * 0.2, 0));
        }
        else
        {
            testValueUpper = (int) round(expectedVal - range * 0.2);
            testValueLower = (int) round(max(expectedVal - range * 1.2, 0));
        }

        int currentValue = (int)((testValueLower + testValueUpper) * 0.5);

        double currentProb = poisson.cumulativeProbability(currentValue);
        double lastProb = currentProb;
        int lastValue = currentValue;
        double probDiff = 0;

        while(iterations < maxIterations)
        {
            probDiff = abs(requiredProb - currentProb) / requiredProb;

            if(probDiff < 0.001)
                return currentValue;

            // if prob is too high, need to lower the test value
            if(currentProb > requiredProb)
            {
                if(currentValue <= testValueLower + 1)
                    break;

                testValueUpper = currentValue;
                currentValue = (int)round((currentValue + testValueLower) * 0.5);
            }
            else
            {
                if(currentValue >= testValueUpper - 1)
                    break;

                testValueLower = currentValue;
                currentValue = (int)round((currentValue + testValueUpper) * 0.5);
            }

            currentProb = poisson.cumulativeProbability(currentValue);
            ++iterations;
        }

        if(iterations >= maxIterations)
        {
            NE_LOGGER.warn(format("max iterations reached: value(%d) test(%d) prob(%.4f diff=%.4f)",
                    expectedVal, currentValue, currentProb, probDiff));
        }

        if(lastProb > requiredProb && currentProb < requiredProb)
        {
            double lastValueFraction = (requiredProb - currentProb) / (lastProb - currentProb);
            return (1 - lastValueFraction) * currentValue + lastValueFraction * lastValue;
        }
        else if(lastProb < requiredProb && currentProb > requiredProb)
        {
            // last prob lower than the current one
            double currentValueFraction = (requiredProb - lastProb) / (currentProb - lastProb);
            return currentValueFraction * currentValue + (1 - currentValueFraction) * lastValue;
        }
        else
        {
            return currentValue;
        }
    }

}
