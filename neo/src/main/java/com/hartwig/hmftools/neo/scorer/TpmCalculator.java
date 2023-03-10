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
import static com.hartwig.hmftools.neo.scorer.NeoRnaData.NO_TPM_VALUE;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
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

    public static final double LOW_PROBABILITY = 0.25;
    public static final double HIGH_PROBABILITY = 1 - LOW_PROBABILITY;

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

        NE_LOGGER.debug("sample({}) TPM norm-factor({})", sampleId, format("%.3f", tpmNormalisationFactor));

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

        // find shared up transcripts across neopitopes from the same variant and cache their total down TPMs
        int totalPeptideCount = 0;

        for(List<NeoEpitopeData> variantNeos : variantNeoMap.values())
        {
            totalPeptideCount += processVariantNeos(variantNeos, tpmNormalisationFactor);
        }

        NE_LOGGER.debug("sample({}) neoepitopes({}) derived {} peptides", sampleId, neoDataList.size(), totalPeptideCount);
    }

    private class SharedUpstreamGroup
    {
        public final Set<Integer> NeIds;
        public final double TpmUp;
        public double TpmDownTotal;

        public SharedUpstreamGroup(final int neId, final double tpmUp, final double tpmDown)
        {
            NeIds = Sets.newHashSet(neId);
            TpmUp = tpmUp;
            TpmDownTotal = tpmDown;
        }

        public String toString() { return format("neIds(%s) tpmUp(%4.3e) tpmTotalDown(%4.3e)", NeIds, TpmUp, TpmDownTotal); }
    }

    private int processVariantNeos(final List<NeoEpitopeData> variantNeos, final double tpmNormalisationFactor)
    {
        // look for neos with shared upstream transcripts
        List<SharedUpstreamGroup> sharedUpstreamGroups = null;

        if(variantNeos.size() > 1)
        {
            for(int i = 0; i < variantNeos.size() - 1; ++i)
            {
                NeoEpitopeData neoData1 = variantNeos.get(i);

                if(sharedUpstreamGroups != null && sharedUpstreamGroups.stream().anyMatch(x -> x.NeIds.contains(neoData1)))
                    continue;

                SharedUpstreamGroup group = null;

                for(int j = i + 1; j < variantNeos.size(); ++j)
                {
                    NeoEpitopeData neoData2 = variantNeos.get(j);

                    if(neoData1.matchingTranscripts(neoData2, FS_UP))
                    {
                        if(group == null)
                        {
                            group = new SharedUpstreamGroup(neoData1.Id, neoData1.RnaData.getTPM(FS_UP), neoData1.RnaData.getTPM(FS_DOWN));

                            if(sharedUpstreamGroups == null)
                                sharedUpstreamGroups = Lists.newArrayList(group);
                            else
                                sharedUpstreamGroups.add(group);
                        }

                        group.NeIds.add(neoData2.Id);
                        group.TpmDownTotal += neoData2.RnaData.getTPM(FS_DOWN);
                    }
                }
            }
        }

        int peptideCount = 0;

        // evaluate peptides per variants to find repeated peptides from different transcripts
        List<PeptideScoreData> peptideScoreDataList = Lists.newArrayList();

        for(NeoEpitopeData neoData : variantNeos)
        {
            List<PeptideData> peptides = EpitopeUtils.generatePeptides(
                    neoData.UpAminoAcids, neoData.NovelAminoAcids, neoData.DownAminoAcids, mPeptideLengthRange, mFlankCount);

            for(PeptideData peptide : peptides)
            {
                if(peptide.Peptide.contains(AA_SELENOCYSTEINE))
                    continue;

                PeptideScoreData peptideScoreData = peptideScoreDataList.stream()
                        .filter(x -> x.Peptide.equals(peptide.Peptide)).findFirst().orElse(null);

                if(peptideScoreData != null)
                    peptideScoreData.addNeoId(neoData.Id);
                else
                    peptideScoreDataList.add(new PeptideScoreData(peptide, neoData.Id));
            }
        }

        peptideCount += peptideScoreDataList.size();

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
                // skip if already calculated
                if(peptideData.tpmCalculated())
                    continue;

                // calculate the TPM Up as total TPM Ups across those neos without double-counting shared upstream transcript neos
                double tpmUp = 0;

                if(peptideData.neIds().size() == 1 && sharedUpstreamGroups == null)
                {
                    tpmUp = neoData.RnaData.getTPM(FS_UP);
                }
                else
                {
                    Set<SharedUpstreamGroup> processedGroups = Sets.newHashSet();

                    for(Integer neId : peptideData.neIds())
                    {
                        SharedUpstreamGroup group = sharedUpstreamGroups != null ?
                                sharedUpstreamGroups.stream().filter(x -> x.NeIds.contains(neId)).findFirst().orElse(null) : null;

                        if(group != null)
                        {
                            if(!processedGroups.contains(group))
                            {
                                processedGroups.add(group);

                                // allocate TPM up based on which shared neos this peptide belongs to
                                double tpmDown = 0;
                                for(int groupNeId : group.NeIds)
                                {
                                    if(!peptideData.neIds().contains(groupNeId))
                                        continue;

                                    NeoEpitopeData peptideNeo = variantNeos.stream().filter(x -> x.Id == groupNeId).findFirst().orElse(null);
                                    tpmDown += peptideNeo.RnaData.getTPM(FS_DOWN);
                                }

                                tpmUp += group.TpmUp * tpmDown / group.TpmDownTotal;
                            }
                        }
                        else
                        {
                            NeoEpitopeData peptideNeo = variantNeos.stream().filter(x -> x.Id == neId).findFirst().orElse(null);
                            tpmUp += peptideNeo.RnaData.getTPM(FS_UP);
                        }
                    }
                }

                // EffectiveTPM = effectiveTPM * [1 – subClonalLikelihood*(1-(min(1,variantCN)]
                double scLikelihood = neoData.VariantType.isPointMutation() ? neoData.SubclonalLikelihood : 1; // note use for fusions
                double cnFactor = 1 - scLikelihood * (1 - min(1.0, neoData.CopyNumber));
                double expectedTpm = tpmUp * cnFactor;

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

        return peptideCount;
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

        if(fragmentCount == 0 && requiredProb == HIGH_PROBABILITY)
            return 0;

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
        else if(fragmentCount == 0)
        {
            currentValue = 1;
            testValueUpper = currentValue * 2;
            testValueLower = fragmentCount;
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
