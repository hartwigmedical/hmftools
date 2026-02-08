package com.hartwig.hmftools.neo.score;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.codon.AminoAcidRna.AA_SELENOCYSTEINE;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.MIN_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.BindConstants.REF_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.score.NeoRnaData.NO_TPM_VALUE;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.neo.NeoEpitopeType;
import com.hartwig.hmftools.common.stats.PoissonCalcs;
import com.hartwig.hmftools.neo.PeptideData;
import com.hartwig.hmftools.neo.epitope.EpitopeUtils;

public class TpmCalculator
{
    private final int mFlankCount;
    private final int[] mPeptideLengthRange;
    private final Map<Integer,PoissonRangeValues> mPoissonRangeValues;

    public static final double LOW_PROBABILITY = 0.25;
    public static final double HIGH_PROBABILITY = 1 - LOW_PROBABILITY;
    public static final double LOW_PROBABILITY_ZERO_OBS_MEAN = 1.39;

    public static final double EFFECTIVE_TPM_ACTUAL_PERC = 0.8;

    public static final int[] DEFAULT_PEPTIDE_LENGTH_RANGE = new int[] { MIN_PEPTIDE_LENGTH, REF_PEPTIDE_LENGTH };

    public TpmCalculator(final int flankCount, final int[] peptideLengthRange)
    {
        mFlankCount = flankCount;
        mPeptideLengthRange = peptideLengthRange;
        mPoissonRangeValues = Maps.newHashMap();
        populatePoissonCache();
    }

    public void compute(final String sampleId, final List<NeoEpitopeData> neoDataList, double samplePloidy)
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
            totalPeptideCount += processVariantNeos(variantNeos, tpmNormalisationFactor, samplePloidy);
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

    private int processVariantNeos(final List<NeoEpitopeData> variantNeos, double tpmNormalisationFactor, double samplePloidy)
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

            PoissonRangeValues rangeValues = getOrCalcRangeValues(neoData.RnaData.fragmentSupport());

            neoData.setCalculatedTpmValues(tpmNormalisationFactor, samplePloidy, rangeValues.LowValue, rangeValues.HighValue);

            double rawEffectiveTpm = neoData.rawEffectiveTpm();

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

                                if(group.TpmDownTotal > 0)
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
                double expectedTpm = tpmUp * neoData.copyNumberFactor();

                double effectiveTpm;

                if(tpmNormalisationFactor == NO_TPM_VALUE || !neoData.RnaData.hasExpression())
                {
                    effectiveTpm = expectedTpm;
                }
                else
                {
                    effectiveTpm = calcEffectiveTpm(
                            tpmNormalisationFactor, rawEffectiveTpm, expectedTpm, rangeValues.LowValue, rangeValues.HighValue);
                }

                peptideData.setEffectiveTpm(effectiveTpm);
            }
        }

        return peptideCount;
    }

    public static double calcEffectiveTpm(
            final double tpmNormalisationFactor, final double rawEffectiveTpm, final double expectedTpm,
            final double probLowValue, final double probHighValue)
    {
        // effectiveTPM = 0.8 * rawEffectiveTPM + 0.2 * max[poisson5%(rawEffectiveTPM), min[Poisson95%(rawEffectiveTPM), expectedTPM]]
        double rangeLow = probLowValue * tpmNormalisationFactor;
        double rangeHigh = probHighValue * tpmNormalisationFactor;

        return EFFECTIVE_TPM_ACTUAL_PERC * rawEffectiveTpm
                + (1 - EFFECTIVE_TPM_ACTUAL_PERC) * max(rangeHigh, min(rangeLow, expectedTpm));
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
            String[] values = line.split(CSV_DELIM, 3);

            mPoissonRangeValues.put(
                    Integer.parseInt( values[0]), new PoissonRangeValues(Double.parseDouble(values[1]), Double.parseDouble(values[2])));
        }
    }

    public static double calcPoissonMean(int fragmentCount, double requiredProb)
    {
        if(fragmentCount < 0)
            return 0;

        if(fragmentCount == 0)
        {
            return requiredProb == HIGH_PROBABILITY ? 0 : LOW_PROBABILITY_ZERO_OBS_MEAN;
        }

        return PoissonCalcs.calcPoissonNoiseValue(fragmentCount, requiredProb);
    }
}
