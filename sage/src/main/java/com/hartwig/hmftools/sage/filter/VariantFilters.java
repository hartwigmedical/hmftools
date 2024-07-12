package com.hartwig.hmftools.sage.filter;

import static java.lang.Math.log10;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;

import static com.hartwig.hmftools.sage.SageConstants.HOTSPOT_MIN_TUMOR_ALT_SUPPORT_SKIP_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.HOTSPOT_MIN_TUMOR_VAF_SKIP_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.HOTSPOT_MIN_ALT_BASE_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.LONG_GERMLINE_INSERT_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MAX_INDEL_GERMLINE_ALT_SUPPORT;
import static com.hartwig.hmftools.sage.SageConstants.MAX_MAP_QUAL_ALT_VS_REF;
import static com.hartwig.hmftools.sage.SageConstants.MAX_READ_EDGE_DISTANCE_PERC;
import static com.hartwig.hmftools.sage.SageConstants.MAX_READ_EDGE_DISTANCE_PROB;
import static com.hartwig.hmftools.sage.SageConstants.STRAND_BIAS_CHECK_THRESHOLD;
import static com.hartwig.hmftools.sage.SageConstants.STRAND_BIAS_NON_ALT_MIN_BIAS;
import static com.hartwig.hmftools.sage.SageConstants.VAF_PROBABILITY_THRESHOLD;
import static com.hartwig.hmftools.sage.SageConstants.VAF_PROBABILITY_THRESHOLD_HOTSPOT;
import static com.hartwig.hmftools.sage.SageConstants.STRAND_BIAS_NON_ALT_MIN_DEPTH;
import static com.hartwig.hmftools.sage.SageConstants.MAP_QUAL_FACTOR_FIXED_PENALTY;
import static com.hartwig.hmftools.sage.filter.SoftFilterConfig.getTieredSoftFilterConfig;

import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.qual.BaseQualAdjustment;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.VariantTier;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class VariantFilters
{
    private final FilterConfig mConfig;
    private final int mReadEdgeDistanceThreshold;

    private final int[] mFilterCounts;

    private enum HardFilterType
    {
        RAW_ALT_SUPPORT,
        TUMOR_QUAL,
        TUMOR_VAF,
        JITTER;
    }

    private static final StrandBiasCalcs mStrandBiasCalcs = new StrandBiasCalcs();

    public VariantFilters(final SageConfig config)
    {
        mConfig = config.Filter;
        mReadEdgeDistanceThreshold = (int)(config.getReadLength() * MAX_READ_EDGE_DISTANCE_PERC);
        mFilterCounts = new int[HardFilterType.values().length];
    }

    public boolean passesHardFilters(final ReadContextCounter readCounter)
    {
        if(readCounter.tier().equals(VariantTier.HOTSPOT))
            return true;

        if(readCounter.altSupport() < mConfig.HardMinTumorRawAltSupport)
        {
            ++mFilterCounts[HardFilterType.RAW_ALT_SUPPORT.ordinal()];
            return false;
        }

        if(readCounter.tumorQuality() < mConfig.HardMinTumorQual)
        {
            ++mFilterCounts[HardFilterType.TUMOR_QUAL.ordinal()];
            return false;
        }

        if(readCounter.vaf() < mConfig.HardMinTumorVaf)
        {
            ++mFilterCounts[HardFilterType.TUMOR_VAF.ordinal()];
            return false;
        }

        if(readCounter.jitter().hardFilterOnNoise())
        {
            ++mFilterCounts[HardFilterType.JITTER.ordinal()];
            return false;
        }

        return true;
    }

    public String filterCountsStr()
    {
        return String.format("alt=%d qual=%d vaf=%d jitter=%d",
                mFilterCounts[HardFilterType.RAW_ALT_SUPPORT.ordinal()],
                mFilterCounts[HardFilterType.TUMOR_QUAL.ordinal()],
                mFilterCounts[HardFilterType.TUMOR_VAF.ordinal()],
                mFilterCounts[HardFilterType.JITTER.ordinal()]);
    }

    public boolean enabled() { return !mConfig.DisableSoftFilter; }

    public void applySoftFilters(final SageVariant variant)
    {
        if(!enabled())
            return;

        final VariantTier tier = variant.tier();
        final SoftFilterConfig softFilterConfig = getTieredSoftFilterConfig(tier, mConfig);

        final Set<SoftFilter> variantFilters = variant.filters();

        final List<ReadContextCounter> refReadCounters = variant.referenceReadCounters();

        // setting ref sample count to zero disables the tumor-germline filters
        int maxReferenceSamples = min(refReadCounters.size(), mConfig.ReferenceSampleCount);

        // where there are multiple tumor samples, if any of them pass then clear any filters from the others
        for(ReadContextCounter tumorReadContextCounter : variant.tumorReadCounters())
        {
            Set<SoftFilter> tumorFilters = Sets.newHashSet();

            applyTumorFilters(tier, softFilterConfig, tumorReadContextCounter, tumorFilters);

            for(int i = 0; i < maxReferenceSamples; ++i)
            {
                ReadContextCounter referenceCounter = refReadCounters.get(i);
                applyTumorGermlineFilters(tier, softFilterConfig, referenceCounter, tumorReadContextCounter, tumorFilters);
            }

            if(tumorFilters.isEmpty())
            {
                variantFilters.clear();
                break;
            }
            else
            {
                variantFilters.addAll(tumorFilters);
            }
        }
    }

    // tumor-only tests
    public void applyTumorFilters(
            final VariantTier tier, final SoftFilterConfig config, final ReadContextCounter primaryTumor, final Set<SoftFilter> filters)
    {
        if(!skipMinTumorQualTest(tier, primaryTumor))
        {
            if(belowMinTumorQual(config, primaryTumor))
                filters.add(SoftFilter.MIN_TUMOR_QUAL);

            if(belowMinTumorVaf(config, primaryTumor))
                filters.add(SoftFilter.MIN_TUMOR_VAF);
        }

        if(belowMaxEdgeDistance(primaryTumor))
        {
            filters.add(SoftFilter.MAX_EDGE_DISTANCE);
        }

        if(exceedsAltVsRefMapQual(primaryTumor))
        {
            filters.add(SoftFilter.MAP_QUAL_REF_ALT_DIFFERENCE);
        }

        if(tier != VariantTier.HOTSPOT)
        {
            if(primaryTumor.belowMinFragmentCoords())
            {
                filters.add(SoftFilter.FRAGMENT_COORDS);
            }

            if(mStrandBiasCalcs.isDepthBelowProbability(primaryTumor.fragmentStrandBiasAlt(), primaryTumor.fragmentStrandBiasNonAlt()))
            {
                filters.add(SoftFilter.FRAGMENT_STRAND_BIAS);
            }

            if(mStrandBiasCalcs.isDepthBelowProbability(primaryTumor.readStrandBiasAlt(), primaryTumor.readStrandBiasNonAlt())
            || (primaryTumor.isIndel() && mStrandBiasCalcs.allOneSide(primaryTumor.readStrandBiasAlt())))
            {
                filters.add(SoftFilter.READ_STRAND_BIAS);
            }

            if(applyJitterFilter(primaryTumor))
            {
                filters.add(SoftFilter.JITTER);
            }
        }

        if(belowMinAverageBaseQuality(primaryTumor, tier))
        {
            filters.add(SoftFilter.MIN_AVG_BASE_QUALITY);
        }
    }

    private boolean skipMinTumorQualTest(final VariantTier tier, final ReadContextCounter primaryTumor)
    {
        return tier == VariantTier.HOTSPOT
                && primaryTumor.altSupport() >= HOTSPOT_MIN_TUMOR_ALT_SUPPORT_SKIP_QUAL
                && Doubles.greaterOrEqual(primaryTumor.vaf(), HOTSPOT_MIN_TUMOR_VAF_SKIP_QUAL)
                && primaryTumor.altBaseQualityTotal() >= HOTSPOT_MIN_ALT_BASE_QUAL;
    }

    // each of the following filters returns true if a variant does not pass the test
    private static boolean belowMinTumorQual(final SoftFilterConfig config, final ReadContextCounter primaryTumor)
    {
        int depth = primaryTumor.depth();
        int altSupport = primaryTumor.altSupport();
        int strongSupport = primaryTumor.strongAltSupport();

        if(strongSupport == 0)
            return true;

        int qualPerRead = (int)round(primaryTumor.qualCounters().modifiedAltBaseQualityTotal() / strongSupport);
        byte qualPerReadFloored = (byte)max(qualPerRead, 15);
        double readQualProb = BaseQualAdjustment.phredQualToProbability(qualPerReadFloored);

        BinomialDistribution distribution = new BinomialDistribution(depth, readQualProb);

        double prob = 1 - distribution.cumulativeProbability(strongSupport - 1);

        double mapQualFactor = calcMapQualFactor(primaryTumor, depth, altSupport, strongSupport);

        primaryTumor.setTumorQualProbability(prob);
        primaryTumor.setMapQualFactor(mapQualFactor);

        if(mapQualFactor < 0)
            return true;

        double scoreCutoff = primaryTumor.useMsiErrorRate() ? config.QualPScoreIndel : config.QualPScore;
        return prob >= scoreCutoff;
    }

    private static double calcMapQualFactor(final ReadContextCounter primaryTumor, final int depth, final int altSupport, final int strongSupport)
    {
        double avgAltModifiedMapQuality = primaryTumor.qualCounters().altModifiedMapQualityTotal() / (double)strongSupport;

        double avgMapQual = primaryTumor.mapQualityTotal() / depth;
        double avgAltMapQual = primaryTumor.altMapQualityTotal() / altSupport;

        double mapQualDiffPenalty = avgMapQual > avgAltMapQual ? 2 * (avgMapQual - avgAltMapQual) : 0;

        double readStrandBiasPenalty;
        StrandBiasData readStrandBiasAlt = primaryTumor.readStrandBiasAlt();
        StrandBiasData readStrandBiasNonAlt = primaryTumor.readStrandBiasNonAlt();
        if(readStrandBiasAlt.minBias() > STRAND_BIAS_CHECK_THRESHOLD)
        {
            readStrandBiasPenalty = 0;
        }
        else if(readStrandBiasNonAlt.forward() < STRAND_BIAS_NON_ALT_MIN_DEPTH || readStrandBiasNonAlt.reverse() < STRAND_BIAS_NON_ALT_MIN_DEPTH)
        {
            readStrandBiasPenalty = 0;
        }
        else if(readStrandBiasNonAlt.minBias() < STRAND_BIAS_NON_ALT_MIN_BIAS && (readStrandBiasNonAlt.bias() < 0.5) == (readStrandBiasAlt.bias() < 0.5))
        {
            readStrandBiasPenalty = 0;
        }
        else
        {
            BinomialDistribution distribution = new BinomialDistribution(readStrandBiasAlt.depth(), 0.5);
            double probability = 2 * distribution.cumulativeProbability((int)round(readStrandBiasAlt.depth() * readStrandBiasAlt.minBias()));
            readStrandBiasPenalty = -10 * log10(probability);
        }

        double avgEdgeDistance = primaryTumor.readEdgeDistance().avgDistanceFromEdge();
        double avgAltEdgeDistance = primaryTumor.readEdgeDistance().avgAltDistanceFromEdge();
        double edgeDistancePenalty = 0;

        if(avgAltEdgeDistance / avgEdgeDistance < 2 * MAX_READ_EDGE_DISTANCE_PERC)
        {
            edgeDistancePenalty = 10 * altSupport * log10(avgEdgeDistance / max(avgAltEdgeDistance, 1));
        }

        double repeatPenalty = 0;

        if(primaryTumor.readContext().MaxRepeat != null && primaryTumor.readContext().MaxRepeat.repeatLength() > 1
        && primaryTumor.readContext().MaxRepeat.repeatLength() * primaryTumor.readContext().MaxRepeat.Count >= 15)
        {
            repeatPenalty = min(3 * primaryTumor.readContext().MaxRepeat.Count, 24);
        }

        double mapQualFactor = avgAltModifiedMapQuality - MAP_QUAL_FACTOR_FIXED_PENALTY - mapQualDiffPenalty - readStrandBiasPenalty
                - edgeDistancePenalty - repeatPenalty;

        return mapQualFactor;
    }

    private static boolean belowMinTumorVaf(final SoftFilterConfig config, final ReadContextCounter primaryTumor)
    {
        return Doubles.lessThan(primaryTumor.vaf(), config.MinTumorVaf) && belowMinTumorVafProbability(primaryTumor);
    }

    private static boolean belowMinTumorVafProbability(final ReadContextCounter primaryTumor)
    {
        int depth = primaryTumor.depth();

        double baseQualAvg = primaryTumor.baseQualityTotal() / (double)depth;

        int supportCount = primaryTumor.strongAltSupport();

        double altBaseQualAvg = primaryTumor.averageAltBaseQuality();

        // SupportCount * min(AvgBQ[ALT] / AvgBQ[DP], 1)
        int adjustedAltSupportCount = supportCount * (int)round(min(altBaseQualAvg / baseQualAvg, 1));

        // p-score is the alternative=’upper’ pvalue for a binomtest with n=DP, k=adj_alt_supp, p=10^(-AvgBQ[DP] / 10)
        double probability = pow(10, -baseQualAvg/10);

        BinomialDistribution distribution = new BinomialDistribution(depth, probability);

        double pValue = 1.0 - distribution.cumulativeProbability(adjustedAltSupportCount - 1);

        if(primaryTumor.tier() == VariantTier.HOTSPOT)
            return pValue > VAF_PROBABILITY_THRESHOLD_HOTSPOT;
        else
            return pValue > VAF_PROBABILITY_THRESHOLD;
    }

    private boolean belowMinAverageBaseQuality(final ReadContextCounter primaryTumor, final VariantTier tier)
    {
        if(primaryTumor.useMsiErrorRate())
            return false;

        if(tier == VariantTier.HOTSPOT)
            return Doubles.lessThan(primaryTumor.averageAltBaseQuality(), mConfig.MinAvgBaseQualHotspot);
        else
            return Doubles.lessThan(primaryTumor.averageAltBaseQuality(), mConfig.MinAvgBaseQual);
    }

    private boolean applyJitterFilter(final ReadContextCounter primaryTumor)
    {
        if(primaryTumor.readContext().MaxRepeat == null)
            return false;

        return primaryTumor.jitter().filterOnNoise();
    }

    private boolean belowMaxEdgeDistance(final ReadContextCounter primaryTumor)
    {
        if(primaryTumor.isLongInsert())
            return false;

        int altMed = primaryTumor.readEdgeDistance().maxAltDistanceFromEdge();

        if(altMed >= mReadEdgeDistanceThreshold)
            return false;

        int maxMed = primaryTumor.readEdgeDistance().maxDistanceFromEdge();

        // note max MED for all reads * 2 covers scenarios were no reads have the variant centred
        double medProb = pow(2 * altMed / (2.0 * maxMed), primaryTumor.altSupport());

        return medProb < MAX_READ_EDGE_DISTANCE_PROB;
    }

    private boolean exceedsAltVsRefMapQual(final ReadContextCounter primaryTumor)
    {
        double depth = primaryTumor.depth();
        double altSupport = primaryTumor.altSupport();

        if(depth == 0 || altSupport == 0)
            return false;

        double avgMapQuality = primaryTumor.mapQualityTotal() / depth;
        double avgAltMapQuality = primaryTumor.altMapQualityTotal() / altSupport;

        return avgMapQuality - avgAltMapQuality > MAX_MAP_QUAL_ALT_VS_REF;
    }

    // germline and paired tumor-germline tests
    public void applyTumorGermlineFilters(
            final VariantTier tier, final SoftFilterConfig config,
            final ReadContextCounter refCounter, final ReadContextCounter primaryTumor, final Set<SoftFilter> filters)
    {
        if(belowMinGermlineCoverage(config, refCounter))
        {
            filters.add(SoftFilter.MIN_GERMLINE_DEPTH);
        }

        if(aboveMaxGermlineVaf(config, refCounter, primaryTumor))
        {
            filters.add(SoftFilter.MAX_GERMLINE_VAF);
        }

        if(aboveMaxGermlineRelativeQual(config, refCounter, primaryTumor))
        {
            filters.add(SoftFilter.MAX_GERMLINE_RELATIVE_QUAL);
        }

        // MNV Tests
        if(aboveMaxMnvIndelGermlineAltSupport(tier, refCounter))
        {
            filters.add(SoftFilter.MAX_GERMLINE_ALT_SUPPORT);
        }
    }

    private static boolean belowMinGermlineCoverage(final SoftFilterConfig config, final ReadContextCounter refCounter)
    {
        boolean chromosomeIsAllosome = HumanChromosome.contains(refCounter.chromosome())
                && HumanChromosome.fromString(refCounter.chromosome()).isAllosome();

        int minGermlineCoverage = chromosomeIsAllosome ? config.MinGermlineCoverageAllosome : config.MinGermlineCoverage;

        return refCounter.depth() < minGermlineCoverage;
    }

    private static boolean aboveMaxGermlineVaf(
            final SoftFilterConfig config, final ReadContextCounter refCounter, final ReadContextCounter primaryTumor)
    {
        double tumorVaf = primaryTumor.vaf();

        if(tumorVaf == 0)
            return false; // will be handled in tumor filters

        int adjustedRefAltCount = refCounter.readCounts().altSupport() + refCounter.simpleAltMatches();

        if(refCounter.variant().indelLengthAbs() > LONG_GERMLINE_INSERT_LENGTH)
        {
            adjustedRefAltCount += refCounter.jitter().shortened() + refCounter.jitter().lengthened();
        }

        double adjustedRefVaf = adjustedRefAltCount / (double)refCounter.readCounts().Total;
        return Doubles.greaterThan(adjustedRefVaf, config.MaxGermlineVaf);
    }

    private static boolean aboveMaxGermlineRelativeQual(
            final SoftFilterConfig config, final ReadContextCounter refCounter, final ReadContextCounter primaryTumor)
    {
        double tumorQual = primaryTumor.tumorQuality();
        double refQual = refCounter.tumorQuality();

        if(tumorQual == 0)
            return false; // will be handled in tumor filters

        double refTumorQualRatio = refQual / tumorQual;
        return Doubles.greaterThan(refTumorQualRatio, config.MaxGermlineVaf);
    }

    private static boolean aboveMaxMnvIndelGermlineAltSupport(final VariantTier tier, final ReadContextCounter refCounter)
    {
        if(tier == VariantTier.HOTSPOT)
            return false;

        if(refCounter.variant().isMNV() || refCounter.isLongInsert())
        {
            double depth = refCounter.depth();
            double altSupportPerc = depth > 0 ? refCounter.altSupport() / depth : 0;
            return altSupportPerc >= MAX_INDEL_GERMLINE_ALT_SUPPORT;
        }

        return false;
    }

    private static final EnumSet<VariantTier> PANEL_ONLY_TIERS = EnumSet.of(VariantTier.HOTSPOT, VariantTier.PANEL);

    public static boolean checkFinalFilters(
            final SageVariant variant, final Set<Integer> passingPhaseSets, final SageConfig config, boolean panelOnly)
    {
        if(panelOnly && !PANEL_ONLY_TIERS.contains(variant.tier()))
            return false;

        if(variant.isPassing())
            return true;

        if(config.Filter.DisableHardFilter)
            return true;

        if(variant.tier() == VariantTier.HOTSPOT)
            return true;

        // Its not always 100% transparent what's happening with the mixed germline dedup logic unless we keep all the associated records
        if(variant.mixedGermlineImpact() > 0)
            return true;

        if(variant.hasReferenceSamples() && variant.hasTumorSamples() && !MitochondrialChromosome.contains(variant.chromosome())
        && !variant.hasMatchingLps(passingPhaseSets))
        {
            final ReadContextCounter refCounter = variant.referenceReadCounters().get(0);

            if(refCounter.altSupport() > config.Filter.FilteredMaxGermlineAltSupport)
                return false;
        }

        return true;
    }
}
