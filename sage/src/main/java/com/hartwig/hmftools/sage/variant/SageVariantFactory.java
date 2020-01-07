package com.hartwig.hmftools.sage.variant;

import java.util.Collections;
import java.util.List;
import java.util.Set;

import javax.annotation.concurrent.NotThreadSafe;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.FilterConfig;
import com.hartwig.hmftools.sage.config.SoftFilterConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.select.TierSelector;

import org.jetbrains.annotations.NotNull;

@NotThreadSafe
public class SageVariantFactory {

    private final FilterConfig config;
    private final TierSelector tierSelector;

    public SageVariantFactory(@NotNull final FilterConfig config, @NotNull final List<VariantHotspot> hotspots,
            @NotNull final List<GenomeRegion> panelRegions) {
        this.config = config;
        this.tierSelector = new TierSelector(panelRegions, hotspots);
    }

    @NotNull
    public SageVariant create(@NotNull final AltContext normal) {
        final SageVariantTier tier = tierSelector.tier(normal);
        final Set<String> filters = germlineOnlyFilters(normal);

        return new SageVariant(tier, filters, normal, Collections.emptyList());
    }

    @NotNull
    public SageVariant create(@NotNull final AltContext normal, @NotNull final List<AltContext> tumorAltContexts) {

        final SageVariantTier tier = tierSelector.tier(normal);
        final SoftFilterConfig softConfig = config.softConfig(tier);
        final Set<String> filters = filters(tier, softConfig, normal, tumorAltContexts.get(0));

        return new SageVariant(tier, filters, normal, tumorAltContexts);
    }


    @NotNull
    private Set<String> germlineOnlyFilters(@NotNull final AltContext germline) {
        final Set<String> result = Sets.newHashSet();

        if (Doubles.lessThan(germline.primaryReadContext().vaf(), config.minGermlineVaf())) {
            result.add(FilterConfig.MIN_GERMLINE_VAF);
        }

        return result;
    }

    @NotNull
    private Set<String> filters(@NotNull final SageVariantTier tier, @NotNull final SoftFilterConfig config,
            @NotNull final AltContext normal, @NotNull final AltContext primaryTumor) {
        Set<String> result = Sets.newHashSet();
        result.addAll(this.config.tumorFilters(tier, primaryTumor));

        final ReadContextCounter normalCounter = normal.primaryReadContext();
        Chromosome contextChromosome = HumanChromosome.fromString(normal.chromosome());
        int minGermlineCoverage =
                contextChromosome.isAllosome() ? config.minGermlineReadContextCoverageAllosome() : config.minGermlineReadContextCoverage();
        if (normal.primaryReadContext().coverage() < minGermlineCoverage) {
            result.add(SoftFilterConfig.MIN_GERMLINE_DEPTH);
        }

        if (Doubles.greaterThan(normalCounter.vaf(), config.maxGermlineVaf())) {
            result.add(SoftFilterConfig.MAX_GERMLINE_VAF);
        }

        double tumorReadContextSupport = primaryTumor.primaryReadContext().altSupport();
        double germlineReadContextSupport = normal.primaryReadContext().altSupport();
        if (Doubles.positive(tumorReadContextSupport)) {
            if (Doubles.greaterThan(germlineReadContextSupport / tumorReadContextSupport, config.maxGermlineRelativeReadContextCount())) {
                result.add(SoftFilterConfig.MAX_GERMLINE_REL_RCC);
            }
        }

        double tumorQual = primaryTumor.primaryReadContext().tumorQuality();
        double germlineQual = normal.primaryReadContext().tumorQuality();
        if (Doubles.positive(tumorQual)) {
            if (Doubles.greaterThan(germlineQual / tumorQual, config.maxGermlineRelativeQual())) {
                result.add(SoftFilterConfig.MAX_GERMLINE_REL_QUAL);
            }
        }

        return result;
    }
}
