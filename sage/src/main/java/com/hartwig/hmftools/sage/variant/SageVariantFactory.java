package com.hartwig.hmftools.sage.variant;

import java.util.List;
import java.util.Set;

import javax.annotation.concurrent.NotThreadSafe;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.FilterConfig;
import com.hartwig.hmftools.sage.config.SoftFilterConfig;
import com.hartwig.hmftools.sage.context.AltContext;
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
    public SageVariant create(@NotNull final AltContext normal, @NotNull final List<AltContext> tumorAltContexts) {

        final SageVariantTier tier = tierSelector.tier(normal);
        final SoftFilterConfig softConfig = softConfig(tier);
        final Set<String> filters =
                tumorAltContexts.isEmpty() ? Sets.newHashSet() : filters(tier, softConfig, normal, tumorAltContexts.get(0));

        return new SageVariant(tier, filters, normal, tumorAltContexts);
    }

    @NotNull
    private SoftFilterConfig softConfig(@NotNull final SageVariantTier tier) {
        switch (tier) {
            case HOTSPOT:
                return config.softHotspotFilter();
            case PANEL:
                return config.softPanelFilter();
            default:
                return config.softWideFilter();
        }
    }

    @NotNull
    private Set<String> filters(@NotNull final SageVariantTier tier, @NotNull final SoftFilterConfig config,
            @NotNull final AltContext normal, @NotNull final AltContext primaryTumor) {
        Set<String> result = Sets.newHashSet();

        if (!skipMinTumorQualTest(tier, primaryTumor) && primaryTumor.primaryReadContext().quality() < config.minTumorQual()) {
            result.add(SoftFilterConfig.MIN_TUMOR_QUAL);
        }

        if (Doubles.lessThan(primaryTumor.primaryReadContext().readContextVaf(), config.minTumorVaf())) {
            result.add(SoftFilterConfig.MIN_TUMOR_VAF);
        }

        if (normal.primaryReadContext().coverage() < config.minGermlineReadContextCoverage()) {
            result.add(SoftFilterConfig.MIN_GERMLINE_DEPTH);
        }

        if (Doubles.greaterThan(normal.supportVaf(), config.maxGermlineVaf())) {
            result.add(SoftFilterConfig.MAX_GERMLINE_VAF);
        }

        double tumorReadContextSupport = primaryTumor.primaryReadContext().support();
        double germlineReadContextSupport = normal.primaryReadContext().support();
        if (Doubles.positive(tumorReadContextSupport)) {
            if (Doubles.greaterThan(germlineReadContextSupport / tumorReadContextSupport, config.maxGermlineRelativeReadContextCount())) {
                result.add(SoftFilterConfig.MAX_GERMLINE_REL_RCC);
            }
        }

        double tumorQual = primaryTumor.primaryReadContext().quality();
        double germlineQual = normal.primaryReadContext().quality();
        if (Doubles.positive(tumorQual)) {
            if (Doubles.greaterThan(germlineQual / tumorQual, config.maxGermlineRelativeQual())) {
                result.add(SoftFilterConfig.MAX_GERMLINE_REL_QUAL);
            }
        }

        return result;
    }

    private boolean skipMinTumorQualTest(@NotNull final SageVariantTier tier, @NotNull final AltContext primaryTumor) {
        return tier.equals(SageVariantTier.HOTSPOT)
                && primaryTumor.primaryReadContext().support() >= config.hotspotMinTumorReadContextSupportToSkipQualCheck();
    }

}
