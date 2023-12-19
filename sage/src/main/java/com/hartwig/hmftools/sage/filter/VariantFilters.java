package com.hartwig.hmftools.sage.filter;

import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;

import static com.hartwig.hmftools.sage.SageConstants.HOTSPOT_MIN_TUMOR_ALT_SUPPORT_SKIP_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.HOTSPOT_MIN_TUMOR_VAF_SKIP_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.HOTSPOT_MIN_RAW_ALT_BASE_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.JITTER_INDEL_MAX_REPEATS;
import static com.hartwig.hmftools.sage.SageConstants.JITTER_INDEL_VAF_THRESHOLD;
import static com.hartwig.hmftools.sage.SageConstants.JITTER_NON_INDEL_MAX_REPEATS;
import static com.hartwig.hmftools.sage.SageConstants.JITTER_NON_INDEL_VAF_THRESHOLD;
import static com.hartwig.hmftools.sage.SageConstants.LONG_GERMLINE_INSERT_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MAX_INDEL_GERMLINE_ALT_SUPPORT;
import static com.hartwig.hmftools.sage.SageConstants.MAX_READ_EDGE_DISTANCE;
import static com.hartwig.hmftools.sage.SageConstants.MAX_READ_EDGE_DISTANCE_PROB;
import static com.hartwig.hmftools.sage.SageConstants.NORMAL_RAW_ALT_BQ_MAX;
import static com.hartwig.hmftools.sage.SageConstants.VAF_PROBABILITY_THRESHOLD;
import static com.hartwig.hmftools.sage.SageConstants.VAF_PROBABILITY_THRESHOLD_HOTSPOT;

import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.VariantTier;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class VariantFilters
{
    private final FilterConfig mConfig;
    private final boolean mHighDepthMode;

    private final int[] mFilterCounts;

    public static int HARD_FC_RAW_BASE_QUAL = 0;
    public static int HARD_FC_RAW_ALT_SUPPORT = 1;
    public static int HARD_FC_TUMOR_QUAL = 2;
    public static int HARD_FC_TUMOR_VAF = 3;

    private static final StrandBiasCalcs mStrandBiasCalcs = new StrandBiasCalcs();

    public VariantFilters(final SageConfig config)
    {
        mConfig = config.Filter;
        mHighDepthMode = config.Quality.HighBaseMode;
        mFilterCounts = new int[HARD_FC_TUMOR_VAF+1];
    }

    public boolean passesHardFilters(final ReadContextCounter readCounter)
    {
        if(readCounter.tier().equals(VariantTier.HOTSPOT))
            return true;

        if(readCounter.rawAltBaseQualityTotal() < mConfig.HardMinTumorRawBaseQuality)
        {
            ++mFilterCounts[HARD_FC_RAW_BASE_QUAL];
            return false;
        }

        if(readCounter.rawAltSupport() < mConfig.HardMinTumorRawAltSupport)
        {
            ++mFilterCounts[HARD_FC_RAW_ALT_SUPPORT];
            return false;
        }

        if(readCounter.tumorQuality() < mConfig.HardMinTumorQual)
        {
            ++mFilterCounts[HARD_FC_TUMOR_QUAL];
            return false;
        }

        if(readCounter.vaf() < mConfig.HardMinTumorVaf)
        {
            ++mFilterCounts[HARD_FC_TUMOR_VAF];
            return false;
        }

        return true;
    }

    public String filterCountsStr()
    {
        return String.format("bq=%d alt=%d qual=%d vaf=%d",
                mFilterCounts[HARD_FC_RAW_BASE_QUAL], mFilterCounts[HARD_FC_RAW_ALT_SUPPORT],
                mFilterCounts[HARD_FC_TUMOR_QUAL], mFilterCounts[HARD_FC_TUMOR_VAF]);
    }

    public boolean enabled() { return !mConfig.DisableSoftFilter; }

    public void applySoftFilters(final SageVariant variant)
    {
        if(!enabled())
            return;

        final VariantTier tier = variant.tier();
        final SoftFilterConfig softFilterConfig = getTieredSoftFilterConfig(tier);

        final Set<String> variantFilters = variant.filters();

        final List<ReadContextCounter> normalReadCounters = variant.normalReadCounters();

        // setting ref sample count to zero disables the tumor-normal filters
        int maxNormalSamples = min(normalReadCounters.size(), mConfig.ReferenceSampleCount);

        // where there are multiple tumor samples, if any of them pass then clear any filters from the others
        for(ReadContextCounter tumorReadContextCounter : variant.tumorReadCounters())
        {
            final Set<String> tumorFilters = Sets.newHashSet();

            applyTumorFilters(tier, softFilterConfig, tumorReadContextCounter, tumorFilters);

            for(int i = 0; i < maxNormalSamples; ++i)
            {
                ReadContextCounter normal = normalReadCounters.get(i);
                applyTumorNormalFilters(tier, softFilterConfig, normal, tumorReadContextCounter, tumorFilters);
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

    public SoftFilterConfig getTieredSoftFilterConfig(final VariantTier tier)
    {
        switch(tier)
        {
            case HOTSPOT:
                return mConfig.SoftHotspotFilter;
            case PANEL:
                return mConfig.SoftPanelFilter;
            case HIGH_CONFIDENCE:
                return mConfig.SoftHighConfidenceFilter;
            default:
                return mConfig.SoftLowConfidenceFilter;
        }
    }

    // tumor-only tests
    public void applyTumorFilters(
            final VariantTier tier, final SoftFilterConfig config, final ReadContextCounter primaryTumor, final Set<String> filters)
    {
        if(!skipMinTumorQualTest(tier, primaryTumor))
        {
            if(belowMinTumorQual(config, primaryTumor))
                filters.add(SoftFilter.MIN_TUMOR_QUAL.filterName());

            if(belowMinTumorVaf(config, primaryTumor))
                filters.add(SoftFilter.MIN_TUMOR_VAF.filterName());
        }

        if(belowMaxEdgeDistance(primaryTumor))
        {
            filters.add(SoftFilter.MAX_EDGE_DISTANCE.filterName());
        }

        if(tier != VariantTier.HOTSPOT)
        {
            if(primaryTumor.belowMinFragmentCoords())
            {
                filters.add(SoftFilter.FRAGMENT_COORDS.filterName());
            }

            boolean checkRefBias = tier == VariantTier.PANEL;

            if(mStrandBiasCalcs.isDepthBelowProbability(primaryTumor.fragmentStrandBiasAlt(), primaryTumor.fragmentStrandBiasRef(), checkRefBias)
            || mStrandBiasCalcs.isDepthBelowProbability(primaryTumor.readStrandBiasAlt(), primaryTumor.readStrandBiasRef(), checkRefBias))
            {
                filters.add(SoftFilter.STRAND_BIAS.filterName());
            }

            if(failsJitterFilter(primaryTumor))
            {
                filters.add(SoftFilter.JITTER.filterName());
            }
        }

        if(belowMinAverageBaseQuality(primaryTumor, tier))
        {
            filters.add(SoftFilter.MIN_AVG_BASE_QUALITY.filterName());
        }
    }

    private boolean skipMinTumorQualTest(final VariantTier tier, final ReadContextCounter primaryTumor)
    {
        return tier == VariantTier.HOTSPOT
                && primaryTumor.altSupport() >= HOTSPOT_MIN_TUMOR_ALT_SUPPORT_SKIP_QUAL
                && Doubles.greaterOrEqual(primaryTumor.vaf(), HOTSPOT_MIN_TUMOR_VAF_SKIP_QUAL)
                && primaryTumor.rawAltBaseQualityTotal() >= HOTSPOT_MIN_RAW_ALT_BASE_QUAL;
    }

    // each of the following filters returns true if a variant does not pass the test
    private static boolean belowMinTumorQual(final SoftFilterConfig config, final ReadContextCounter primaryTumor)
    {
        return primaryTumor.tumorQuality() < config.MinTumorQual;
    }

    private static boolean belowMinTumorVaf(final SoftFilterConfig config, final ReadContextCounter primaryTumor)
    {
        return Doubles.lessThan(primaryTumor.vaf(), config.MinTumorVaf) && belowMinTumorVafProbability(primaryTumor);
    }

    private static boolean belowMinTumorVafProbability(final ReadContextCounter primaryTumor)
    {
        int rawDepth = primaryTumor.rawDepth();

        double rawBaseQualAvg = (primaryTumor.rawRefBaseQualityTotal() + primaryTumor.rawAltBaseQualityTotal()) / (double)rawDepth;
        double bqrBaseQualAvg = primaryTumor.recalibratedBaseQualityTotal() / (double)rawDepth;
        double minBaseQualAvg = min(rawBaseQualAvg, bqrBaseQualAvg);

        int supportCount = primaryTumor.strongAltSupport();

        double altBaseQualAvg = primaryTumor.averageAltBaseQuality();

        // SupportCount * min(AvgBQ[ALT] / AvgBQ[DP], 1)
        int adjustedAltSupportCount = supportCount * (int)round(min(altBaseQualAvg / minBaseQualAvg, 1));

        // p-score is the alternative=’upper’ pvalue for a binomtest with n=DP, k=adj_alt_supp, p=10^(-AvgBQ[DP] / 10)
        double probability = pow(10, -minBaseQualAvg/10);

        BinomialDistribution distribution = new BinomialDistribution(rawDepth, probability);

        double pValue = 1.0 - distribution.cumulativeProbability(adjustedAltSupportCount - 1);

        if(primaryTumor.tier() == VariantTier.HOTSPOT)
            return pValue > VAF_PROBABILITY_THRESHOLD_HOTSPOT;
        else
            return pValue > VAF_PROBABILITY_THRESHOLD;
    }

    private boolean belowMinAverageBaseQuality(final ReadContextCounter primaryTumor, final VariantTier tier)
    {
        if(tier == VariantTier.HOTSPOT)
            return Doubles.lessThan(primaryTumor.averageAltBaseQuality(), mConfig.MinAvgBaseQualHotspot);
        else
            return Doubles.lessThan(primaryTumor.averageAltBaseQuality(), mConfig.MinAvgBaseQual);
    }

    private static boolean isLongInsert(final ReadContextCounter variant)
    {
        return variant.isIndel() && variant.alt().length() > LONG_GERMLINE_INSERT_LENGTH;
    }

    private boolean failsJitterFilter(final ReadContextCounter primaryTumor)
    {
        if(!mHighDepthMode)
            return false;

        int maxRepeats = primaryTumor.readContext().RepeatCount; // repeat count not available yet

        if(primaryTumor.isIndel())
        {
            // INDELs if inserted/deleted bases == RC_MH and VAF < (MAX_REP - 3) * 0.0125
            String indelBases = primaryTumor.variant().isInsert() ?
                    primaryTumor.alt().substring(1) : primaryTumor.ref().substring(1);

            double vafLimit = (maxRepeats - JITTER_INDEL_MAX_REPEATS) * JITTER_INDEL_VAF_THRESHOLD;
            return indelBases.equals(primaryTumor.readContext().Microhomology) && primaryTumor.vaf() < vafLimit;
        }
        else
        {
            // SNVs/MNVs if MAX_REP > 5 and VAF < 0.01
            return maxRepeats > JITTER_NON_INDEL_MAX_REPEATS && primaryTumor.vaf() < JITTER_NON_INDEL_VAF_THRESHOLD;
        }
    }

    private boolean belowMaxEdgeDistance(final ReadContextCounter primaryTumor)
    {
        if(isLongInsert(primaryTumor))
            return false;

        int altMed = primaryTumor.readEdgeDistance().maxAltDistanceFromUnclippedEdge();
        int maxMed = primaryTumor.readEdgeDistance().maxDistanceFromUnclippedEdge();

        if(altMed >= MAX_READ_EDGE_DISTANCE)
            return false;

        // note max MED for all reads * 2 covers scenarios were no reads have the variant centred
        double medProb = pow(2 * altMed / (2.0 * maxMed), primaryTumor.altSupport());

        return medProb < MAX_READ_EDGE_DISTANCE_PROB;
    }

    // normal and paired tumor-normal tests
    public void applyTumorNormalFilters(
            final VariantTier tier, final SoftFilterConfig config,
            final ReadContextCounter normal, final ReadContextCounter primaryTumor, final Set<String> filters)
    {
        if(belowMinGermlineCoverage(config, normal))
        {
            filters.add(SoftFilter.MIN_GERMLINE_DEPTH.filterName());
        }

        if(aboveMaxGermlineVaf(config, normal, primaryTumor))
        {
            filters.add(SoftFilter.MAX_GERMLINE_VAF.filterName());
        }

        // Paired Tests
        if(aboveMaxGermlineQual(config, normal, primaryTumor))
        {
            filters.add(SoftFilter.MAX_GERMLINE_REL_RAW_BASE_QUAL.filterName());
        }

        // MNV Tests
        if(aboveMaxMnvIndelNormalAltSupport(tier, normal))
        {
            filters.add(SoftFilter.MAX_GERMLINE_ALT_SUPPORT.filterName());
        }
    }

    private static boolean belowMinGermlineCoverage(final SoftFilterConfig config, final ReadContextCounter normal)
    {
        boolean chromosomeIsAllosome = HumanChromosome.contains(normal.chromosome())
                && HumanChromosome.fromString(normal.chromosome()).isAllosome();

        boolean isLongInsert = isLongInsert(normal);

        int minGermlineCoverage = chromosomeIsAllosome ?
                (isLongInsert ? config.MinGermlineCoverageAllosomeLongInsert : config.MinGermlineCoverageAllosome)
                : (isLongInsert ? config.MinGermlineCoverageLongInsert : config.MinGermlineCoverage);

        return normal.depth() < minGermlineCoverage;
    }

    private static boolean aboveMaxGermlineVaf(
            final SoftFilterConfig config, final ReadContextCounter normal, final ReadContextCounter primaryTumor)
    {
        double normalVaf = normal.vaf();

        if(!primaryTumor.isIndel() && normal.rawAltBaseQualityTotal() > 0 && normal.rawAltBaseQualityTotal() < NORMAL_RAW_ALT_BQ_MAX
        && normal.rawAltSupport() == normal.altSupport())
        {
            double normalRawBqVaf = normal.rawAltBaseQualityTotal() / (double)(normal.rawAltBaseQualityTotal() + normal.rawRefBaseQualityTotal());
            normalVaf = min(normalVaf, normalRawBqVaf);
        }

        return Doubles.greaterThan(normalVaf, config.MaxGermlineVaf);
    }

    private static boolean aboveMaxGermlineQual(
            final SoftFilterConfig config, final ReadContextCounter normal, final ReadContextCounter primaryTumor)
    {
        boolean isLongInsert = primaryTumor.isIndel() && normal.alt().length() > LONG_GERMLINE_INSERT_LENGTH;
        double tumorQual = isLongInsert ? primaryTumor.tumorQuality() : primaryTumor.rawAltBaseQualityTotal();
        double germlineQual = isLongInsert ? normal.tumorQuality() : normal.rawAltBaseQualityTotal();

        return Doubles.positive(tumorQual) && Doubles.greaterThan(germlineQual / tumorQual, config.MaxGermlineRelativeQual);
    }

    private static boolean aboveMaxMnvIndelNormalAltSupport(final VariantTier tier, final ReadContextCounter normal)
    {
        if(tier == VariantTier.HOTSPOT)
            return false;

        if(normal.variant().isMNV() || (normal.variant().isInsert() && normal.variant().indelLength() >= LONG_GERMLINE_INSERT_LENGTH))
        {
            double depth = normal.depth();
            double altSupportPerc = depth > 0 ? normal.altSupport() / depth : 0;
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

        if(!variant.isNormalEmpty() && !variant.isTumorEmpty() && !MitochondrialChromosome.contains(variant.chromosome())
                && !variant.hasMatchingLps(passingPhaseSets))
        {
            final ReadContextCounter normal = variant.normalReadCounters().get(0);

            if(normal.altSupport() > config.Filter.FilteredMaxNormalAltSupport)
                return false;
        }

        return true;
    }


}
