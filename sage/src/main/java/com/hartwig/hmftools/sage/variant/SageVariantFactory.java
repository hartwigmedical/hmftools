package com.hartwig.hmftools.sage.variant;

import java.util.List;
import java.util.Set;

import javax.annotation.concurrent.NotThreadSafe;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.sage.config.FilterConfig;
import com.hartwig.hmftools.sage.config.SoftFilterConfig;
import com.hartwig.hmftools.sage.context.AltContext;

import org.jetbrains.annotations.NotNull;

@NotThreadSafe
public class SageVariantFactory {

    private final InPanel inPanel;
    private final FilterConfig config;
    private final InHotspot inHotspot;

    public SageVariantFactory(@NotNull final Chromosome chromosome, @NotNull final FilterConfig config, @NotNull final List<VariantHotspot> hotspots,
            @NotNull final List<GenomeRegion> panelRegions) {
        this.config = config;
        this.inPanel = new InPanel(chromosome, panelRegions);
        this.inHotspot = new InHotspot(chromosome, hotspots);
    }

    public SageVariantFactory(@NotNull final FilterConfig config, @NotNull final ListMultimap<Chromosome, VariantHotspot> hotspots,
            @NotNull final ListMultimap<Chromosome, GenomeRegion> panelRegions) {
        this.config = config;
        this.inPanel = new InPanel(panelRegions);
        this.inHotspot = new InHotspot(hotspots);
    }

    @NotNull
    public SageVariant create(@NotNull final AltContext normal, @NotNull final List<AltContext> tumorAltContexts) {

        final SageVariantTier tier = tier(normal);
        final SoftFilterConfig softConfig = softConfig(tier);
        final Set<String> filters = filters(softConfig, normal, tumorAltContexts.get(0));

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
    private SageVariantTier tier(@NotNull final AltContext normal) {
        if (inHotspot.isOnHotspot(normal)) {
            return SageVariantTier.HOTSPOT;
        }
        if (inPanel.inPanel(normal)) {
            return SageVariantTier.PANEL;
        }
        return SageVariantTier.WIDE;
    }

    @NotNull
    private Set<String> filters(@NotNull final SoftFilterConfig config, @NotNull final AltContext normal,
            @NotNull final AltContext primaryTumor) {
        Set<String> result = Sets.newHashSet();

        if (primaryTumor.primaryReadContext().quality() < config.minTumorQual()) {
            result.add(SoftFilterConfig.MIN_TUMOR_QUAL);
        }

        if (Doubles.lessThan(primaryTumor.primaryReadContext().readContextVaf(), config.minTumorVaf())) {
            result.add(SoftFilterConfig.MIN_TUMOR_VAF);
        }

        if (normal.readDepth() < config.minGermlineDepth()) {
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

}
