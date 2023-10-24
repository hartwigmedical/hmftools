package com.hartwig.hmftools.sage.filter;

import static java.lang.Math.min;

import static com.hartwig.hmftools.sage.SageConstants.HOTSPOT_MIN_TUMOR_ALT_SUPPORT_SKIP_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.HOTSPOT_MIN_TUMOR_VAF_SKIP_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.HOTSPOT_MIN_RAW_ALT_BASE_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.LONG_GERMLINE_INSERT_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MAX_INDEL_GERMLINE_ALT_SUPPORT;
import static com.hartwig.hmftools.sage.SageConstants.NORMAL_RAW_ALT_BQ_MAX;

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

public class VariantFilters
{
    private final FilterConfig mConfig;
    private final StrandBiasCalcs mStrandBiasCalcs;

    private final int[] mFilterCounts;

    public static int HARD_FC_RAW_BASE_QUAL = 0;
    public static int HARD_FC_RAW_ALT_SUPPORT = 1;
    public static int HARD_FC_TUMOR_QUAL = 2;
    public static int HARD_FC_TUMOR_VAF = 3;

    public VariantFilters(final FilterConfig config)
    {
        mConfig = config;
        mStrandBiasCalcs = new StrandBiasCalcs();
        mFilterCounts = new int[HARD_FC_TUMOR_VAF+1];
    }

    public boolean passesHardFilters(final ReadContextCounter readCounter)
    {
        if(readCounter.tier().equals(VariantTier.HOTSPOT))
            return true;

        if(readCounter.rawAltBaseQuality() < mConfig.HardMinTumorRawBaseQuality)
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
        return String.format("bq=%d a=%d tq=%d qv=%d",
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

        if(mStrandBiasCalcs.isDepthBelowProbability(primaryTumor.strandBias(), primaryTumor.strandDepth()))
        {
            filters.add(SoftFilter.STRAND_BIAS.filterName());
        }

        if(belowMinAverageBaseQuality(primaryTumor, tier))
        {
            filters.add(SoftFilter.MIN_AVG_BASE_QUALITY.filterName());
        }
    }

    private boolean skipMinTumorQualTest(final VariantTier tier, final ReadContextCounter primaryTumor)
    {
        return tier.equals(VariantTier.HOTSPOT)
                && primaryTumor.altSupport() >= HOTSPOT_MIN_TUMOR_ALT_SUPPORT_SKIP_QUAL
                && Doubles.greaterOrEqual(primaryTumor.vaf(), HOTSPOT_MIN_TUMOR_VAF_SKIP_QUAL)
                && primaryTumor.rawAltBaseQuality() >= HOTSPOT_MIN_RAW_ALT_BASE_QUAL;
    }

    // each of the following filters returns true if a variant does not pass the test
    private static boolean belowMinTumorQual(final SoftFilterConfig config, final ReadContextCounter primaryTumor)
    {
        return primaryTumor.tumorQuality() < config.MinTumorQual;
    }

    private static boolean belowMinTumorVaf(final SoftFilterConfig config, final ReadContextCounter primaryTumor)
    {
        return Doubles.lessThan(primaryTumor.vaf(), config.MinTumorVaf);
    }

    private boolean belowMinAverageBaseQuality(final ReadContextCounter primaryTumor, final VariantTier tier)
    {
        if(tier == VariantTier.HOTSPOT)
            return Doubles.lessThan(primaryTumor.averageAltBaseQuality(), mConfig.MinAvgBaseQualHotspot);
        else
            return Doubles.lessThan(primaryTumor.averageAltBaseQuality(), mConfig.MinAvgBaseQual);
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

        boolean isLongInsert = normal.isIndel() && normal.alt().length() > LONG_GERMLINE_INSERT_LENGTH;

        int minGermlineCoverage = chromosomeIsAllosome ?
                (isLongInsert ? config.MinGermlineCoverageAllosomeLongInsert : config.MinGermlineCoverageAllosome)
                : (isLongInsert ? config.MinGermlineCoverageLongInsert : config.MinGermlineCoverage);

        return normal.depth() < minGermlineCoverage;
    }

    private static boolean aboveMaxGermlineVaf(
            final SoftFilterConfig config, final ReadContextCounter normal, final ReadContextCounter primaryTumor)
    {
        double normalVaf = normal.vaf();

        if(!primaryTumor.isIndel() && normal.rawAltBaseQuality() > 0 && normal.rawAltBaseQuality() < NORMAL_RAW_ALT_BQ_MAX
        && normal.rawAltSupport() == normal.altSupport())
        {
            double normalRawBqVaf = normal.rawAltBaseQuality() / (double)(normal.rawAltBaseQuality() + normal.rawRefBaseQuality());
            normalVaf = min(normalVaf, normalRawBqVaf);
        }

        return Doubles.greaterThan(normalVaf, config.MaxGermlineVaf);
    }

    private static boolean aboveMaxGermlineQual(
            final SoftFilterConfig config, final ReadContextCounter normal, final ReadContextCounter primaryTumor)
    {
        boolean isLongInsert = primaryTumor.isIndel() && normal.alt().length() > LONG_GERMLINE_INSERT_LENGTH;
        double tumorQual = isLongInsert ? primaryTumor.tumorQuality() : primaryTumor.rawAltBaseQuality();
        double germlineQual = isLongInsert ? normal.tumorQuality() : normal.rawAltBaseQuality();

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
