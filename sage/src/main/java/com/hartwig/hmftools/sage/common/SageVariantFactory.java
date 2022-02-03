package com.hartwig.hmftools.sage.common;

import static java.lang.Math.min;

import static com.hartwig.hmftools.sage.SageConstants.HOTSPOT_MIN_TUMOR_ALT_SUPPORT_SKIP_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.HOTSPOT_MIN_TUMOR_VAF_SKIP_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.NORMAL_RAW_ALT_BQ_MAX;

import java.util.List;
import java.util.Set;

import javax.annotation.concurrent.NotThreadSafe;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.config.FilterConfig;
import com.hartwig.hmftools.sage.config.SoftFilter;
import com.hartwig.hmftools.sage.config.SoftFilterConfig;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

@NotThreadSafe
public class SageVariantFactory
{
    private final FilterConfig mConfig;

    public SageVariantFactory(final FilterConfig config)
    {
        mConfig = config;
    }

    public SageVariant create(final Candidate candidate, final List<ReadContextCounter> normal, final List<ReadContextCounter> tumor)
    {
        boolean isNormalEmpty = normal.isEmpty();

        final VariantTier tier = candidate.tier();
        final SoftFilterConfig softConfig = mConfig.softConfig(tier);

        SageVariant variant = new SageVariant(candidate, normal, tumor);

        if(!mConfig.SoftFilter)
        {
            return new SageVariant(candidate, normal, tumor);
        }

        final Set<String> variantFilters = variant.filters();

        for(ReadContextCounter tumorReadContextCounter : tumor)
        {
            final Set<String> tumorFilters = Sets.newHashSet();

            tumorFilters(tier, softConfig, tumorReadContextCounter, tumorFilters);

            if(!isNormalEmpty)
            {
                somaticFilters(tier, softConfig, normal.get(0), tumorReadContextCounter, tumorFilters);
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

        return new SageVariant(candidate, normal, tumor);
    }

    private void tumorFilters(
            final VariantTier tier, final SoftFilterConfig config, final ReadContextCounter primaryTumor, final Set<String> filters)
    {
        // TUMOR Tests
        final boolean skipTumorTests = skipMinTumorQualTest(tier, primaryTumor);
        if(!skipTumorTests && primaryTumor.tumorQuality() < config.MinTumorQual)
        {
            filters.add(SoftFilter.MIN_TUMOR_QUAL.toString());
        }

        if(!skipTumorTests && Doubles.lessThan(primaryTumor.vaf(), config.MinTumorVaf))
        {
            filters.add(SoftFilter.MIN_TUMOR_VAF.toString());
        }
    }

    private void somaticFilters(
            final VariantTier tier, final SoftFilterConfig config,
            final ReadContextCounter normal, final ReadContextCounter primaryTumor, final Set<String> filters)
    {
        // Germline Tests
        boolean chromosomeIsAllosome = HumanChromosome.contains(normal.chromosome())
                && HumanChromosome.fromString(normal.chromosome()).isAllosome();

        int minGermlineCoverage = chromosomeIsAllosome ? config.MinGermlineReadContextCoverageAllosome : config.MinGermlineReadContextCoverage;

        if(normal.coverage() < minGermlineCoverage)
        {
            filters.add(SoftFilter.MIN_GERMLINE_DEPTH.toString());
        }

        double normalVaf = normal.vaf();

        if(!primaryTumor.isIndel() && normal.rawAltBaseQuality() > 0 && normal.rawAltBaseQuality() < NORMAL_RAW_ALT_BQ_MAX)
        {
            double normalRawBqVcf = normal.rawAltBaseQuality() / (double)(normal.rawAltBaseQuality() + normal.rawRefBaseQuality());
            normalVaf = min(normalVaf, normalRawBqVcf);
        }

        if(Doubles.greaterThan(normalVaf, config.MaxGermlineVaf))
        {
            filters.add(SoftFilter.MAX_GERMLINE_VAF.toString());
        }

        // Paired Tests
        double tumorQual = primaryTumor.rawAltBaseQuality();
        double germlineQual = normal.rawAltBaseQuality();
        if(Doubles.positive(tumorQual))
        {
            if(Doubles.greaterThan(germlineQual / tumorQual, config.MaxGermlineRelativeQual))
            {
                filters.add(SoftFilter.MAX_GERMLINE_REL_RAW_BASE_QUAL.toString());
            }
        }

        // MNV Tests
        if(tier != VariantTier.HOTSPOT && normal.variant().isMNV() && this.mConfig.MnvFilter)
        {
            if(normal.altSupport() != 0)
            {
                filters.add(SoftFilter.MAX_GERMLINE_ALT_SUPPORT.toString());
            }
        }
    }

    private boolean skipMinTumorQualTest(final VariantTier tier, final ReadContextCounter primaryTumor)
    {
        return tier.equals(VariantTier.HOTSPOT) && primaryTumor.altSupport() >= HOTSPOT_MIN_TUMOR_ALT_SUPPORT_SKIP_QUAL
                && Doubles.greaterOrEqual(primaryTumor.vaf(), HOTSPOT_MIN_TUMOR_VAF_SKIP_QUAL);
    }
}
