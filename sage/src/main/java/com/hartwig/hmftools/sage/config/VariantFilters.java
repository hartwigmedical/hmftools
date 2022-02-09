package com.hartwig.hmftools.sage.config;

import static java.lang.Math.min;

import static com.hartwig.hmftools.sage.SageConstants.HOTSPOT_MIN_TUMOR_ALT_SUPPORT_SKIP_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.HOTSPOT_MIN_TUMOR_VAF_SKIP_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.NORMAL_RAW_ALT_BQ_MAX;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.VariantTier;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

public class VariantFilters
{
    private final FilterConfig mConfig;

    public VariantFilters(final FilterConfig config)
    {
        mConfig = config;
    }

    public boolean passesHardFilters(final ReadContextCounter readContextCounter)
    {
        if(readContextCounter.Tier.equals(VariantTier.HOTSPOT))
            return true;

        return readContextCounter.rawAltBaseQuality() >= mConfig.HardMinTumorRawBaseQuality
                && readContextCounter.rawAltSupport() >= mConfig.HardMinTumorRawAltSupport
                && readContextCounter.tumorQuality() >= mConfig.HardMinTumorQual;
    }

    public boolean enabled() { return mConfig.SoftFilter; }

    public void applySoftFilters(final SageVariant variant)
    {
        if(!enabled())
            return;

        final VariantTier tier = variant.tier();
        final SoftFilterConfig softFilterConfig = getTieredSoftFilterConfig(tier);

        final Set<String> variantFilters = variant.filters();

        final ReadContextCounter normal = !variant.isNormalEmpty() ? variant.normalReadCounters().get(0) : null;

        // where there are multiple tumor samples, if any of them pass then clear any filters from the others
        for(ReadContextCounter tumorReadContextCounter : variant.tumorReadCounters())
        {
            final Set<String> tumorFilters = Sets.newHashSet();

            applyTumorFilters(tier, softFilterConfig, tumorReadContextCounter, tumorFilters);

            if(normal != null)
            {
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
                filters.add(SoftFilter.MIN_TUMOR_QUAL.toString());

            if(belowMinTumorVaf(config, primaryTumor))
                filters.add(SoftFilter.MIN_TUMOR_VAF.toString());
        }
    }

    private boolean skipMinTumorQualTest(final VariantTier tier, final ReadContextCounter primaryTumor)
    {
        return tier.equals(VariantTier.HOTSPOT)
                && primaryTumor.altSupport() >= HOTSPOT_MIN_TUMOR_ALT_SUPPORT_SKIP_QUAL
                && Doubles.greaterOrEqual(primaryTumor.vaf(), HOTSPOT_MIN_TUMOR_VAF_SKIP_QUAL);
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

    // normal and paired tumor-normal tests
    public void applyTumorNormalFilters(
            final VariantTier tier, final SoftFilterConfig config,
            final ReadContextCounter normal, final ReadContextCounter primaryTumor, final Set<String> filters)
    {
        if(belowMinGermlineCoverage(config, normal))
        {
            filters.add(SoftFilter.MIN_GERMLINE_DEPTH.toString());
        }

        if(aboveMaxGermlineVaf(config, normal, primaryTumor))
        {
            filters.add(SoftFilter.MAX_GERMLINE_VAF.toString());
        }

        // Paired Tests
        if(aboveMaxGermlineQual(config, normal, primaryTumor))
        {
            filters.add(SoftFilter.MAX_GERMLINE_REL_RAW_BASE_QUAL.toString());
        }

        // MNV Tests
        if(aboveMaxMnvNormalAltSupport(tier, normal, this.mConfig.MnvFilter))
        {
            filters.add(SoftFilter.MAX_GERMLINE_ALT_SUPPORT.toString());
        }
    }

    private static boolean belowMinGermlineCoverage(final SoftFilterConfig config, final ReadContextCounter normal)
    {
        boolean chromosomeIsAllosome = HumanChromosome.contains(normal.chromosome())
                && HumanChromosome.fromString(normal.chromosome()).isAllosome();

        int minGermlineCoverage = chromosomeIsAllosome ? config.MinGermlineReadContextCoverageAllosome : config.MinGermlineReadContextCoverage;

        return normal.coverage() < minGermlineCoverage;
    }

    private static boolean aboveMaxGermlineVaf(
            final SoftFilterConfig config, final ReadContextCounter normal, final ReadContextCounter primaryTumor)
    {
        double normalVaf = normal.vaf();

        if(!primaryTumor.isIndel() && normal.rawAltBaseQuality() > 0 && normal.rawAltBaseQuality() < NORMAL_RAW_ALT_BQ_MAX)
        {
            double normalRawBqVcf = normal.rawAltBaseQuality() / (double)(normal.rawAltBaseQuality() + normal.rawRefBaseQuality());
            normalVaf = min(normalVaf, normalRawBqVcf);
        }

        return Doubles.greaterThan(normalVaf, config.MaxGermlineVaf);
    }

    private static boolean aboveMaxGermlineQual(
            final SoftFilterConfig config, final ReadContextCounter normal, final ReadContextCounter primaryTumor)
    {
        // Paired Tests
        double tumorQual = primaryTumor.rawAltBaseQuality();
        double germlineQual = normal.rawAltBaseQuality();

        return Doubles.positive(tumorQual) && Doubles.greaterThan(germlineQual / tumorQual, config.MaxGermlineRelativeQual);
    }

    private static boolean aboveMaxMnvNormalAltSupport(
            final VariantTier tier, final ReadContextCounter normal, boolean applyMnvFilter)
    {
        return tier != VariantTier.HOTSPOT && normal.variant().isMNV() && applyMnvFilter && normal.altSupport() != 0;
    }

    /*
    public boolean isSoftFiltered(final ReadContextCounter normal, final ReadContextCounter tumor)
    {
        if(!mConfig.SoftFilter)
            return false;

        SoftFilterConfig softFilterConfig = getTieredSoftFilterConfig(tumor.Tier);

        if(!skipMinTumorQualTest(tumor.Tier, tumor))
        {
            if(belowMinTumorQual(softFilterConfig, tumor))
                return true;

            if(belowMinTumorVaf(softFilterConfig, tumor))
                return true;
        }

        if(belowMinGermlineCoverage(softFilterConfig, normal))
            return true;

        if(aboveMaxGermlineVaf(softFilterConfig, normal, tumor))
            return true;

        // Paired Tests
        if(aboveMaxGermlineQual(softFilterConfig, normal, tumor))
            return true;

        // MNV Tests
        if(aboveMaxMnvNormalAltSupport(tumor.Tier, normal, this.mConfig.MnvFilter))
            return true;

        return false;
    }
    */

}
