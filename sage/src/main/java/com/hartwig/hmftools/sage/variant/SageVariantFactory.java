package com.hartwig.hmftools.sage.variant;

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
import com.hartwig.hmftools.sage.read.ReadContextCounter;

import org.jetbrains.annotations.NotNull;

@NotThreadSafe
public class SageVariantFactory
{
    private final FilterConfig mConfig;

    public SageVariantFactory(@NotNull final FilterConfig config)
    {
        this.mConfig = config;
    }

    @NotNull
    public SageVariant create(@NotNull final Candidate candidate, @NotNull final List<ReadContextCounter> normal,
            @NotNull final List<ReadContextCounter> tumor)
    {
        assert (!tumor.isEmpty());
        final Set<String> filters = Sets.newHashSet();
        final boolean isNormalEmpty = normal.isEmpty();

        final SageVariantTier tier = candidate.tier();
        final SoftFilterConfig softConfig = mConfig.softConfig(tier);

        if(!mConfig.SoftFilter)
        {
            return new SageVariant(candidate, filters, normal, tumor);
        }

        boolean passingTumor = false;
        for(ReadContextCounter tumorReadContextCounter : tumor)
        {
            final Set<String> tumorFilters = Sets.newHashSet();

            tumorFilters.addAll(tumorFilters(tier, softConfig, tumorReadContextCounter));
            if(!isNormalEmpty)
            {
                tumorFilters.addAll(somaticFilters(tier, softConfig, normal.get(0), tumorReadContextCounter));
            }

            if(tumorFilters.isEmpty())
            {
                passingTumor = true;
            }

            filters.addAll(tumorFilters);
        }

        if(passingTumor)
        {
            filters.clear();
        }

        return new SageVariant(candidate, filters, normal, tumor);
    }

    @NotNull
    private Set<String> tumorFilters(@NotNull final SageVariantTier tier, @NotNull final SoftFilterConfig config,
            @NotNull final ReadContextCounter primaryTumor)
    {
        final Set<String> result = Sets.newHashSet();

        // TUMOR Tests
        final boolean skipTumorTests = skipMinTumorQualTest(tier, primaryTumor);
        if(!skipTumorTests && primaryTumor.tumorQuality() < config.MinTumorQual)
        {
            result.add(SoftFilter.MIN_TUMOR_QUAL.toString());
        }

        if(!skipTumorTests && Doubles.lessThan(primaryTumor.vaf(), config.MinTumorVaf))
        {
            result.add(SoftFilter.MIN_TUMOR_VAF.toString());
        }

        return result;
    }

    @NotNull
    private Set<String> somaticFilters(@NotNull final SageVariantTier tier, @NotNull final SoftFilterConfig config,
            @NotNull final ReadContextCounter normal, @NotNull final ReadContextCounter primaryTumor)
    {
        Set<String> result = Sets.newHashSet();

        // Germline Tests
        final boolean chromosomeIsAllosome =
                HumanChromosome.contains(normal.chromosome()) && HumanChromosome.fromString(normal.chromosome()).isAllosome();
        int minGermlineCoverage =
                chromosomeIsAllosome ? config.MinGermlineReadContextCoverageAllosome : config.MinGermlineReadContextCoverage;
        if(normal.coverage() < minGermlineCoverage)
        {
            result.add(SoftFilter.MIN_GERMLINE_DEPTH.toString());
        }

        if(Doubles.greaterThan(normal.vaf(), config.MaxGermlineVaf))
        {
            result.add(SoftFilter.MAX_GERMLINE_VAF.toString());
        }

        // Paired Tests
        double tumorQual = primaryTumor.rawAltBaseQuality();
        double germlineQual = normal.rawAltBaseQuality();
        if(Doubles.positive(tumorQual))
        {
            if(Doubles.greaterThan(germlineQual / tumorQual, config.MaxGermlineRelativeQual))
            {
                result.add(SoftFilter.MAX_GERMLINE_REL_RAW_BASE_QUAL.toString());
            }
        }

        // MNV Tests
        if(tier != SageVariantTier.HOTSPOT && normal.variant().isMNV() && this.mConfig.MnvFilter)
        {
            if(normal.altSupport() != 0)
            {
                result.add(SoftFilter.MAX_GERMLINE_ALT_SUPPORT.toString());
            }
        }

        return result;
    }

    private boolean skipMinTumorQualTest(@NotNull final SageVariantTier tier, @NotNull final ReadContextCounter primaryTumor)
    {
        return tier.equals(SageVariantTier.HOTSPOT) && primaryTumor.altSupport() >= mConfig.hotspotMinTumorAltSupportToSkipQualCheck()
                && Doubles.greaterOrEqual(primaryTumor.vaf(), mConfig.hotspotMinTumorVafToSkipQualCheck());
    }
}
