package com.hartwig.hmftools.sage.filter;

import static java.lang.Math.log10;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.VariantTier.HOTSPOT;
import static com.hartwig.hmftools.common.variant.VariantTier.PANEL;
import static com.hartwig.hmftools.sage.ReferenceData.isHighlyPolymorphic;
import static com.hartwig.hmftools.sage.SageConfig.isSbx;
import static com.hartwig.hmftools.sage.SageConfig.isUltima;
import static com.hartwig.hmftools.sage.SageConfig.isIllumina;
import static com.hartwig.hmftools.sage.SageConstants.HIGHLY_POLYMORPHIC_GENES_ALT_MAP_QUAL_THRESHOLD;
import static com.hartwig.hmftools.sage.SageConstants.HOTSPOT_MIN_TUMOR_ALT_SUPPORT_SKIP_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.HOTSPOT_MIN_TUMOR_VAF_SKIP_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.HOTSPOT_MIN_ALT_BASE_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.MAP_QUAL_INDEL_REPEAT_PENALTY;
import static com.hartwig.hmftools.sage.SageConstants.MAP_QUAL_NON_INDEL_REPEAT_PENALTY;
import static com.hartwig.hmftools.sage.SageConstants.MAP_QUAL_READ_BIAS_CAP;
import static com.hartwig.hmftools.sage.SageConstants.MAX_GERMLINE_VAF_PANEL_INDEL_REPEAT_THRESHOLD_FACTOR;
import static com.hartwig.hmftools.sage.SageConstants.MAX_GERMLINE_VAF_PANEL_INDEL_REPEAT_VAF_FACTOR;
import static com.hartwig.hmftools.sage.SageConstants.MAX_GERMLINE_VAF_PANEL_VAF_FACTOR;
import static com.hartwig.hmftools.sage.SageConstants.MAX_GERMLINE_QUAL_HET_TUMOR_VAF;
import static com.hartwig.hmftools.sage.SageConstants.MAX_GERMLINE_QUAL_PROB_HOTSPOT;
import static com.hartwig.hmftools.sage.SageConstants.MAX_GERMLINE_QUAL_PROB_OTHER;
import static com.hartwig.hmftools.sage.SageConstants.MAX_GERMLINE_QUAL_PROB_PANEL;
import static com.hartwig.hmftools.sage.SageConstants.MAX_GERMLINE_QUAL_RATIO_THRESHOLD;
import static com.hartwig.hmftools.sage.SageConstants.MAX_GERMLINE_QUAL_RATIO_THRESHOLD_HOTSPOT;
import static com.hartwig.hmftools.sage.SageConstants.MAX_GERMLINE_VAF_THRESHOLD_FACTOR;
import static com.hartwig.hmftools.sage.SageConstants.MAX_INDEL_GERMLINE_ALT_SUPPORT;
import static com.hartwig.hmftools.sage.SageConstants.MAX_MAP_QUAL_ALT_VS_REF;
import static com.hartwig.hmftools.sage.SageConstants.MAX_READ_EDGE_DISTANCE_PERC;
import static com.hartwig.hmftools.sage.SageConstants.MAX_READ_EDGE_DISTANCE_PERC_PANEL;
import static com.hartwig.hmftools.sage.SageConstants.MAX_READ_EDGE_DISTANCE_PROB;
import static com.hartwig.hmftools.sage.SageConstants.MIN_TQP_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.MIN_TQP_QUAL_MSI_VARIANT;
import static com.hartwig.hmftools.sage.SageConstants.REALIGNED_MAX_PERC;
import static com.hartwig.hmftools.sage.SageConstants.REQUIRED_STRONG_SUPPORT;
import static com.hartwig.hmftools.sage.SageConstants.REQUIRED_STRONG_SUPPORT_HOTSPOT;
import static com.hartwig.hmftools.sage.SageConstants.REQUIRED_UNIQUE_FRAG_COORDS_1;
import static com.hartwig.hmftools.sage.SageConstants.REQUIRED_UNIQUE_FRAG_COORDS_2;
import static com.hartwig.hmftools.sage.SageConstants.REQUIRED_UNIQUE_FRAG_COORDS_AD_1;
import static com.hartwig.hmftools.sage.SageConstants.REQUIRED_UNIQUE_FRAG_COORDS_AD_2;
import static com.hartwig.hmftools.sage.SageConstants.STRAND_BIAS_CHECK_THRESHOLD;
import static com.hartwig.hmftools.sage.SageConstants.STRAND_BIAS_NON_ALT_MIN_BIAS;
import static com.hartwig.hmftools.sage.SageConstants.VAF_PROBABILITY_THRESHOLD;
import static com.hartwig.hmftools.sage.SageConstants.VAF_PROBABILITY_THRESHOLD_HOTSPOT;
import static com.hartwig.hmftools.sage.SageConstants.STRAND_BIAS_NON_ALT_MIN_DEPTH;
import static com.hartwig.hmftools.sage.SageConstants.MAP_QUAL_FACTOR_FIXED_PENALTY;
import static com.hartwig.hmftools.sage.SageConstants.HIGHLY_POLYMORPHIC_GENES_MAX_QUALITY;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_BASE_QUAL_FIXED_PENALTY;
import static com.hartwig.hmftools.sage.SageConstants.GERMLINE_HET_MIN_EXPECTED_VAF;
import static com.hartwig.hmftools.sage.SageConstants.GERMLINE_HET_MIN_SAMPLING_PROB;
import static com.hartwig.hmftools.sage.filter.SoftFilterConfig.getTieredSoftFilterConfig;
import static com.hartwig.hmftools.sage.seqtech.SbxUtils.MQF_NM_1_THRESHOLD_DEDUCTION;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.belowExpectedHpQuals;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.belowExpectedT0Quals;

import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.redux.BaseQualAdjustment;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class VariantFilters
{
    private final boolean mIsGermline;
    private final FilterConfig mConfig;

    private final int[] mFilterCounts;

    private enum HardFilterType
    {
        RAW_ALT_SUPPORT,
        TUMOR_QUAL,
        TUMOR_VAF,
        JITTER;
    }

    public static final StrandBiasCalcs STRAND_BIAS_CALCS = new StrandBiasCalcs();

    public VariantFilters(final SageConfig config)
    {
        mConfig = config.Filter;
        mIsGermline = config.IsGermline;
        mFilterCounts = new int[HardFilterType.values().length];
    }

    public static double readEdgeDistanceThreshold(final VariantTier tier)
    {
        return tier == HOTSPOT || tier == PANEL ? MAX_READ_EDGE_DISTANCE_PERC_PANEL : MAX_READ_EDGE_DISTANCE_PERC;
    }

    public boolean passesHardFilters(final ReadContextCounter readCounter)
    {
        if(readCounter.tier().equals(HOTSPOT))
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

        VariantTier tier = variant.tier();
        SoftFilterConfig softFilterConfig = getTieredSoftFilterConfig(tier, mConfig);

        Set<SoftFilter> variantFilters = variant.filters();

        List<ReadContextCounter> refReadCounters = variant.referenceReadCounters();

        // setting ref sample count to zero disables the tumor-germline filters
        int maxReferenceSamples = min(refReadCounters.size(), mConfig.ReferenceSampleCount);

        // where there are multiple tumor samples, if any of them pass then clear any filters from the others
        for(ReadContextCounter tumorReadContextCounter : variant.tumorReadCounters())
        {
            Set<SoftFilter> tumorFilters = Sets.newHashSet();
            applyTumorFilters(variant, tier, softFilterConfig, tumorReadContextCounter, tumorFilters);

            if(!mIsGermline)
            {
                for(int i = 0; i < maxReferenceSamples; ++i)
                {
                    ReadContextCounter referenceCounter = refReadCounters.get(i);
                    applyTumorGermlineFilters(tier, softFilterConfig, referenceCounter, tumorReadContextCounter, tumorFilters);
                }
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
            final SageVariant variant, final VariantTier tier, final SoftFilterConfig config, final ReadContextCounter primaryTumor,
            final Set<SoftFilter> filters)
    {
        if(!skipMinTumorQualTest(tier, primaryTumor))
        {
            if(belowMinTumorQual(config, tier, primaryTumor, mIsGermline))
                filters.add(SoftFilter.MIN_TUMOR_QUAL);

            if(belowMinMapQualFactor(variant, config, tier, primaryTumor))
                filters.add(SoftFilter.MIN_MAP_QUAL_FACTOR);

            if(belowMinTumorVaf(config, primaryTumor))
                filters.add(SoftFilter.MIN_TUMOR_VAF);
        }

        if(belowMaxEdgeDistance(tier, primaryTumor))
        {
            filters.add(SoftFilter.MAX_EDGE_DISTANCE);
        }

        if(exceedsAltVsRefMapQual(primaryTumor))
        {
            filters.add(SoftFilter.MAP_QUAL_REF_ALT_DIFFERENCE);
        }

        if(belowMinFragmentCoords(primaryTumor))
        {
            filters.add(SoftFilter.FRAGMENT_COORDS);
        }

        if(belowMinStrongSupport(primaryTumor))
        {
            filters.add(SoftFilter.MIN_TUMOR_SUPPORT);
        }

        /* NOTE: disabled due to impact on ctDNA samples
        if(exceedsAltFragmentLength(primaryTumor))
        {
            filters.add(SoftFilter.FRAGMENT_LENGTH);
        }
        */

        if(exceedsRealignedPercentage(primaryTumor))
        {
            filters.add(SoftFilter.REALIGNED_FREQ);
        }

        if(tier != HOTSPOT)
        {
            if(STRAND_BIAS_CALCS.isDepthBelowProbability(primaryTumor.fragmentStrandBiasAlt(), primaryTumor.fragmentStrandBiasNonAlt()))
            {
                filters.add(SoftFilter.FRAGMENT_STRAND_BIAS);
            }

            if(STRAND_BIAS_CALCS.isDepthBelowProbability(primaryTumor.readStrandBiasAlt(), primaryTumor.readStrandBiasNonAlt())
            || (primaryTumor.isIndel() && STRAND_BIAS_CALCS.allOneSide(primaryTumor.readStrandBiasAlt())))
            {
                filters.add(SoftFilter.READ_STRAND_BIAS);
            }
        }

        if(applyJitterFilter(primaryTumor))
        {
            filters.add(SoftFilter.JITTER);
        }

        if(isUltima())
        {
            if(belowExpectedHpQuals(primaryTumor))
            {
                filters.add(SoftFilter.MIN_AVG_HP_QUAL);
            }

            if(belowExpectedT0Quals(primaryTumor))
            {
                filters.add(SoftFilter.MIN_AVG_T0_QUAL);
            }
        }
        else if(belowMinAverageBaseQuality(primaryTumor, tier))
        {
            filters.add(SoftFilter.MIN_AVG_BASE_QUALITY);
        }
    }

    private boolean skipMinTumorQualTest(final VariantTier tier, final ReadContextCounter primaryTumor)
    {
        return tier == HOTSPOT
                && primaryTumor.altSupport() >= HOTSPOT_MIN_TUMOR_ALT_SUPPORT_SKIP_QUAL
                && Doubles.greaterOrEqual(primaryTumor.vaf(), HOTSPOT_MIN_TUMOR_VAF_SKIP_QUAL)
                && primaryTumor.qualCounters().altRecalibratedBaseQualityTotal() >= HOTSPOT_MIN_ALT_BASE_QUAL;
    }

    private static boolean boostNovelIndel(final VariantTier tier, final ReadContextCounter primaryTumor)
    {
        return (tier == HOTSPOT || tier == PANEL || isIllumina())
                && primaryTumor.isIndel() && !primaryTumor.useMsiErrorRate()
                && primaryTumor.jitter().lengthened() + primaryTumor.jitter().shortened() == 0
                && (double)primaryTumor.nonAltNmCountTotal() / (primaryTumor.depth() - primaryTumor.altSupport()) < 0.75;
    }

    // each of the following filters returns true if a variant does not pass the test
    private static boolean belowMinTumorQual(
            final SoftFilterConfig config, final VariantTier tier, final ReadContextCounter primaryTumor, final boolean isGermline)
    {
        int depth = primaryTumor.depth();
        int altSupport = primaryTumor.altSupport();
        int strongSupport = primaryTumor.strongAltSupport();

        if(strongSupport == 0)
        {
            primaryTumor.setTumorQualProbability(1.0);
            return true;
        }

        int strongNonMediumSupport = strongSupport - primaryTumor.strongMediumQualSupport();
        double recalibratedAltMediumBaseQualityTotal = primaryTumor.qualCounters().strongAltMediumRecalibratedBaseQualityTotal();
        double recalibratedAltNonMediumBaseQualityTotal = primaryTumor.qualCounters().strongAltRecalibratedBaseQualityTotal() - recalibratedAltMediumBaseQualityTotal;
        int mediumSupportContribution = (int)(strongNonMediumSupport * recalibratedAltMediumBaseQualityTotal/recalibratedAltNonMediumBaseQualityTotal);
        int adjustedStrongSupport = strongNonMediumSupport + mediumSupportContribution;
        double altNonMediumFinalBaseQualityTotal = primaryTumor.qualCounters().altFinalBaseQualityTotal() - primaryTumor.qualCounters().altMediumFinalBaseQualityTotal();

        int qualPerRead = (int)round(altNonMediumFinalBaseQualityTotal / strongNonMediumSupport);

        if(boostNovelIndel(tier, primaryTumor))
            qualPerRead += DEFAULT_BASE_QUAL_FIXED_PENALTY;  // should boost by the actual config base qual penalty

        int minQualForTqp = primaryTumor.qualCache().isMsiSampleAndVariant() ? MIN_TQP_QUAL_MSI_VARIANT : MIN_TQP_QUAL;

        byte qualPerReadFloored = (byte)max(qualPerRead, minQualForTqp);
        double readQualProb = BaseQualAdjustment.phredQualToProbability(qualPerReadFloored);

        BinomialDistribution distribution = new BinomialDistribution(depth, readQualProb);

        double prob = 1 - distribution.cumulativeProbability(adjustedStrongSupport - 1);

        if(isGermline)
        {
            BinomialDistribution hetGermlineDistribution = new BinomialDistribution(depth, GERMLINE_HET_MIN_EXPECTED_VAF);
            double hetProb = hetGermlineDistribution.cumulativeProbability(altSupport);
            if(hetProb < prob)
            {
                prob /= max(hetProb, GERMLINE_HET_MIN_SAMPLING_PROB);
                prob = min(prob, 1.0);
            }
        }

        primaryTumor.setTumorQualProbability(prob);

        double scoreCutoff = config.QualPScore;
        return prob >= scoreCutoff;
    }

    private static boolean belowMinMapQualFactor(
            final SageVariant variant, final SoftFilterConfig config, final VariantTier tier, final ReadContextCounter primaryTumor)
    {
        int depth = primaryTumor.depth();
        int altSupport = primaryTumor.altSupport();

        if(depth == 0 || altSupport == 0)
            return false;

        double mapQualFactor = calcMapQualFactor(variant, tier, primaryTumor, depth, altSupport, primaryTumor.strongAltSupport());

        primaryTumor.setMapQualFactor(mapQualFactor);

        double threshold = config.MapQualFactor;

        if(isSbx() && primaryTumor.minNumberOfEvents() == 1)
            threshold -= MQF_NM_1_THRESHOLD_DEDUCTION;

        return mapQualFactor < threshold;
    }

    private static double calcMapQualFactor(
            final SageVariant variant, final VariantTier tier, final ReadContextCounter primaryTumor,
            int depth, int altSupport, int strongSupport)
    {
        double avgAltFinalMapQuality = primaryTumor.qualCounters().altFinalMapQualityTotal() / (double)strongSupport;

        double avgMapQual = primaryTumor.mapQualityTotal() / (double)depth;
        double avgAltMapQual = primaryTumor.altMapQualityTotal() / (double)altSupport;

        double mapQualDiffPenalty = max(2 * (avgMapQual - avgAltMapQual), 0);

        boolean highlyPolymorphicSite = isHighlyPolymorphic(primaryTumor.variant());

        if(highlyPolymorphicSite)
        {
            avgAltFinalMapQuality = min(MAP_QUAL_FACTOR_FIXED_PENALTY + HIGHLY_POLYMORPHIC_GENES_MAX_QUALITY, avgAltMapQual);

            if(avgAltMapQual >= HIGHLY_POLYMORPHIC_GENES_ALT_MAP_QUAL_THRESHOLD)
                mapQualDiffPenalty = 0;
        }

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

            if(probability > 0)
                readStrandBiasPenalty = -10 * log10(probability);
            else
                readStrandBiasPenalty = MAP_QUAL_READ_BIAS_CAP; // fall-back to apply a penalty for extreme bias
        }

        double avgEdgeDistance = primaryTumor.readEdgeDistance().avgDistanceFromEdge();
        double avgAltEdgeDistance = primaryTumor.readEdgeDistance().avgAltDistanceFromEdge();
        double edgeDistancePenalty = 0;

        double readEdgeDistanceThresholdPerc = readEdgeDistanceThreshold(tier);
        double altAvgEdgeDistanceRatio = avgAltEdgeDistance / avgEdgeDistance;

        if(!highlyPolymorphicSite && altAvgEdgeDistanceRatio < 2 * readEdgeDistanceThresholdPerc && !primaryTumor.isLongIndel())
        {
            edgeDistancePenalty = 10 * altSupport * log10(avgEdgeDistance / max(avgAltEdgeDistance, 0.001));

            if(!isSbx() && !variant.nearMultiBaseIndel())
                edgeDistancePenalty = min(edgeDistancePenalty, 10);
        }

        double repeatPenalty = 0;

        if(primaryTumor.readContext().MaxRepeat != null && primaryTumor.readContext().MaxRepeat.repeatLength() > 1
        && primaryTumor.readContext().MaxRepeat.repeatLength() * primaryTumor.readContext().MaxRepeat.Count >= 15)
        {
            int maxPenalty = primaryTumor.isIndel() ? MAP_QUAL_INDEL_REPEAT_PENALTY : MAP_QUAL_NON_INDEL_REPEAT_PENALTY;
            repeatPenalty = min(3 * primaryTumor.readContext().MaxRepeat.Count, maxPenalty);
        }

        double mapQualFactor = avgAltFinalMapQuality - MAP_QUAL_FACTOR_FIXED_PENALTY - mapQualDiffPenalty - readStrandBiasPenalty
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

        double baseQualAvg = primaryTumor.qualCounters().recalibratedBaseQualityTotal() / (double)depth;

        int supportCount = primaryTumor.strongAltSupport();

        double altBaseQualAvg = primaryTumor.averageAltRecalibratedBaseQuality();

        if(boostNovelIndel(primaryTumor.tier(), primaryTumor))
            altBaseQualAvg += DEFAULT_BASE_QUAL_FIXED_PENALTY;

        // SupportCount * min(AvgBQ[ALT] / AvgBQ[DP], 1)
        int adjustedAltSupportCount = supportCount * (int)round(min(altBaseQualAvg / baseQualAvg, 1));

        // p-score is the alternative=’upper’ pvalue for a binomtest with n=DP, k=adj_alt_supp, p=10^(-AvgBQ[DP] / 10)
        double probability = pow(10, -baseQualAvg/10);

        BinomialDistribution distribution = new BinomialDistribution(depth, probability);

        double pValue = 1.0 - distribution.cumulativeProbability(adjustedAltSupportCount - 1);

        if(primaryTumor.tier() == HOTSPOT)
            return pValue > VAF_PROBABILITY_THRESHOLD_HOTSPOT;
        else
            return pValue > VAF_PROBABILITY_THRESHOLD;
    }

    private boolean belowMinAverageBaseQuality(final ReadContextCounter primaryTumor, final VariantTier tier)
    {
        int threshold = tier == HOTSPOT ? mConfig.MinAvgBaseQualHotspot : mConfig.MinAvgBaseQual;

        if(primaryTumor.useMsiErrorRate())
            return false;

        double avgBaseQuality = primaryTumor.averageAltSeqTechBaseQuality();

        return Doubles.lessThan(avgBaseQuality, threshold);
    }

    private boolean applyJitterFilter(final ReadContextCounter primaryTumor)
    {
        if(primaryTumor.readContext().MaxRepeat == null)
            return false;

        return primaryTumor.jitter().filterOnNoise();
    }

    private boolean belowMaxEdgeDistance(final VariantTier tier, final ReadContextCounter primaryTumor)
    {
        if(primaryTumor.isLongInsert())
            return false;

        double altMed = primaryTumor.readEdgeDistance().maxAltDistanceFromEdge();

        double edgeDistanceThreshold = readEdgeDistanceThreshold(tier);

        if(altMed >= edgeDistanceThreshold)
            return false;

        double maxMed = primaryTumor.readEdgeDistance().maxDistanceFromEdge();

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

    private boolean belowMinStrongSupport(final ReadContextCounter primaryTumor)
    {
        int strongSupportThreshold = primaryTumor.tier() == HOTSPOT ? REQUIRED_STRONG_SUPPORT_HOTSPOT : REQUIRED_STRONG_SUPPORT;
        int strongSupport;
        if(primaryTumor.useMsiErrorRate())
            strongSupport = primaryTumor.jitter().validQualFullSupport();
        else if(isUltima())
            strongSupport = primaryTumor.strongAltSupport();  // depending on HP lengths, may not have any high qual reads
        else
            strongSupport = primaryTumor.strongHighQualSupport();
        return strongSupport < strongSupportThreshold;
    }

    private boolean belowMinFragmentCoords(final ReadContextCounter primaryTumor)
    {
        if(primaryTumor.fragmentCoords() == null)
            return false;

        int minRequiredUniqueFrags = 1;
        if(primaryTumor.altSupport() >= REQUIRED_UNIQUE_FRAG_COORDS_AD_2)
            minRequiredUniqueFrags = REQUIRED_UNIQUE_FRAG_COORDS_2;
        else if(primaryTumor.altSupport() >= REQUIRED_UNIQUE_FRAG_COORDS_AD_1)
            minRequiredUniqueFrags = REQUIRED_UNIQUE_FRAG_COORDS_1;

        return primaryTumor.fragmentCoords().minCount() < minRequiredUniqueFrags;
    }

    private boolean exceedsRealignedPercentage(final ReadContextCounter primaryTumor)
    {
        if(!primaryTumor.isLongIndel())
        {
            int altSupport = primaryTumor.altSupport();

            if(altSupport == 0)
                return false;

            double realignedPerc = primaryTumor.readCounts().Realigned / (double)altSupport;

            if(realignedPerc > REALIGNED_MAX_PERC)
                return true;
        }

        return false;
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

        if(aboveMaxGermlineVaf(tier, config, refCounter, primaryTumor))
        {
            filters.add(SoftFilter.MAX_GERMLINE_VAF);
        }

        if(aboveMaxGermlineRelativeQual(tier, refCounter, primaryTumor))
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
            final VariantTier tier, final SoftFilterConfig config, final ReadContextCounter refCounter,
            final ReadContextCounter primaryTumor)
    {
        double tumorVaf = primaryTumor.vaf();

        if(tumorVaf == 0)
            return false; // will be handled in tumor filters

        int adjustedRefAltCount = refCounter.readCounts().altSupport() + refCounter.simpleAltMatches();

        if(refCounter.isLongIndel())
        {
            adjustedRefAltCount += refCounter.jitter().shortened() + refCounter.jitter().lengthened();
        }

        boolean isPanelIndelRepeat = isPanelIndelRepeatVariant(primaryTumor);
        return aboveMaxGermlineVaf(tier, isPanelIndelRepeat, tumorVaf, adjustedRefAltCount, refCounter.readCounts().Total, config.MaxGermlineVaf);
    }

    public static boolean isPanelIndelRepeatVariant(final ReadContextCounter counter)
    {
        return counter.tier() == PANEL && counter.variant().isIndel() && counter.readContext().MaxRepeat != null;
    }

    public static boolean aboveMaxGermlineVaf(
            final VariantTier tier, boolean isPanelIndelRepeat, double tumorVaf, int adjustedRefAltCount, int refTotalReads,
            double maxGermlineVaf)
    {
        if(tumorVaf == 0)
            return false; // will be handled in tumor filters

        if(tier == PANEL || tier == HOTSPOT)
        {
            double threshold = tumorVaf;

            if(tier == PANEL)
            {
                if(isPanelIndelRepeat)
                    threshold /= MAX_GERMLINE_VAF_PANEL_INDEL_REPEAT_VAF_FACTOR;
                else
                    threshold /= MAX_GERMLINE_VAF_PANEL_VAF_FACTOR;
            }

            double minThresholdFactor = isPanelIndelRepeat ?
                    MAX_GERMLINE_VAF_PANEL_INDEL_REPEAT_THRESHOLD_FACTOR : MAX_GERMLINE_VAF_THRESHOLD_FACTOR;

            maxGermlineVaf = max(min(threshold, maxGermlineVaf * minThresholdFactor), maxGermlineVaf);
        }

        double adjustedRefVaf = adjustedRefAltCount / (double)refTotalReads;
        return Doubles.greaterThan(adjustedRefVaf, maxGermlineVaf);
    }

    private static boolean aboveMaxGermlineRelativeQual(
            final VariantTier tier, final ReadContextCounter refCounter, final ReadContextCounter primaryTumor)
    {
        int tumorAvgAltBaseQuality = (int)round(primaryTumor.averageAltRecalibratedBaseQuality());
        int refAvgAltBaseQuality = (int)round(refCounter.averageAltRecalibratedBaseQuality());
        int tumorQuality = primaryTumor.readQuals().Full + primaryTumor.readQuals().PartialCore + primaryTumor.readQuals().Realigned;
        int refQuality = refCounter.readQuals().Full + refCounter.readQuals().PartialCore + refCounter.readQuals().Realigned;

        return aboveMaxGermlineRelativeQual(
                tier, tumorQuality, primaryTumor.vaf(), primaryTumor.depth(), tumorAvgAltBaseQuality,
                refQuality, refCounter.depth(), refCounter.altSupport(), refAvgAltBaseQuality);
    }

    public static boolean aboveMaxGermlineRelativeQual(
            final VariantTier tier,
            double tumorQual, double tumorVaf, int tumorDepth, double tumorAvgBaseQual,
            double refQual, int refDepth, int refAltSupport, double refAvgBaseQual)
    {
        if(tumorQual == 0)
            return false; // will be handled in tumor filters

        if(refDepth == 0 || tumorDepth == 0)
            return false; // rely on other germline filters, eg min germline depth

        // adjust inputs by avg base qual
        double adjustedTumorVaf = min(max(tumorVaf - pow(10, -tumorAvgBaseQual / 10), 0.0), MAX_GERMLINE_QUAL_HET_TUMOR_VAF);

        int adjustedRefAd = max((int)round(refAltSupport - refDepth * pow(10, -refAvgBaseQual / 10)), 0);
        double depthRatio = refDepth / (double)tumorDepth;
        double adjustedRefQualRatio = refQual / tumorQual / depthRatio;

        BinomialDistribution distribution = new BinomialDistribution(refDepth, adjustedTumorVaf);

        double prob = distribution.cumulativeProbability(adjustedRefAd);

        double probThreshold = tier == HOTSPOT ? MAX_GERMLINE_QUAL_PROB_HOTSPOT :
                (tier == PANEL ? MAX_GERMLINE_QUAL_PROB_PANEL : MAX_GERMLINE_QUAL_PROB_OTHER);

        double ratioThreshold = tier == HOTSPOT ? MAX_GERMLINE_QUAL_RATIO_THRESHOLD_HOTSPOT : MAX_GERMLINE_QUAL_RATIO_THRESHOLD;

        return prob > probThreshold && adjustedRefQualRatio > ratioThreshold;
    }

    private static boolean aboveMaxMnvIndelGermlineAltSupport(final VariantTier tier, final ReadContextCounter refCounter)
    {
        return aboveMaxMnvIndelGermlineAltSupport(
                tier, refCounter.variant().isMNV(), refCounter.isLongInsert(), refCounter.depth(), refCounter.altSupport());
    }

    public static boolean aboveMaxMnvIndelGermlineAltSupport(
            final VariantTier tier, boolean isMnv, boolean isLongInsert, int refDepth, int refAltSupport)
    {
        if(tier == HOTSPOT)
            return false;

        if(isMnv || isLongInsert)
        {
            double altSupportPerc = refDepth > 0 ? refAltSupport / (double)refDepth : 0;
            return altSupportPerc >= MAX_INDEL_GERMLINE_ALT_SUPPORT;
        }

        return false;
    }

    private static final EnumSet<VariantTier> PANEL_ONLY_TIERS = EnumSet.of(HOTSPOT, PANEL);

    public static boolean checkFinalFilters(
            final SageVariant variant, final Set<Integer> passingPhaseSets, final SageConfig config, boolean panelOnly)
    {
        if(panelOnly && !PANEL_ONLY_TIERS.contains(variant.tier()))
            return false;

        if(variant.isPassing())
            return true;

        if(config.Filter.DisableHardFilter)
            return true;

        if(variant.tier() == HOTSPOT || variant.tier() == PANEL)
            return true;

        // Its not always 100% transparent what's happening with the mixed germline dedup logic unless we keep all the associated records
        if(variant.mixedGermlineImpact() > 0)
            return true;

        if(!config.IsGermline
        && variant.hasReferenceSamples() && variant.hasTumorSamples() && !MitochondrialChromosome.contains(variant.chromosome())
        && !variant.hasMatchingLps(passingPhaseSets))
        {
            final ReadContextCounter refCounter = variant.referenceReadCounters().get(0);

            if(refCounter.altSupport() > config.Filter.FilteredMaxGermlineAltSupport)
                return false;
        }

        return true;
    }
}
