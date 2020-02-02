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
import com.hartwig.hmftools.sage.config.SoftFilter;
import com.hartwig.hmftools.sage.config.SoftFilterConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.select.TierSelector;
import com.hartwig.hmftools.sage.vcf.SageVCF;

import org.jetbrains.annotations.NotNull;

@NotThreadSafe
public class SageVariantFactory {

    private final FilterConfig config;
    private final TierSelector tierSelector;

    public SageVariantFactory(@NotNull final FilterConfig config, @NotNull final List<VariantHotspot> hotspots,
            @NotNull final List<GenomeRegion> panelRegions, @NotNull final List<GenomeRegion> highConfidenceRegions) {
        this.config = config;
        this.tierSelector = new TierSelector(hotspots, panelRegions, highConfidenceRegions);
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
        final Set<String> filters = pairedFilters(tier, softConfig, normal, tumorAltContexts.get(0));

        return new SageVariant(tier, filters, normal, tumorAltContexts);
    }

    @NotNull
    private Set<String> germlineOnlyFilters(@NotNull final AltContext germline) {
        final Set<String> result = Sets.newHashSet();

        if (Doubles.lessThan(germline.primaryReadContext().vaf(), config.minGermlineVaf())) {
            result.add(SageVCF.MIN_GERMLINE_VAF);
        }

        return result;
    }

    @NotNull
    private Set<String> pairedFilters(@NotNull final SageVariantTier tier, @NotNull final SoftFilterConfig config,
            @NotNull final AltContext normal, @NotNull final AltContext primaryTumor) {
        Set<String> result = Sets.newHashSet();

        // TUMOR Tests
        final boolean skipTumorTests = skipMinTumorQualTest(tier, primaryTumor);
        if (!skipTumorTests && primaryTumor.primaryReadContext().tumorQuality() < config.minTumorQual()) {
            result.add(SoftFilter.MIN_TUMOR_QUAL.toString());
        }

        if (!skipTumorTests && Doubles.lessThan(primaryTumor.primaryReadContext().vaf(), config.minTumorVaf())) {
            result.add(SoftFilter.MIN_TUMOR_VAF.toString());
        }

        // GERMLINE Tests
        final ReadContextCounter normalCounter = normal.primaryReadContext();
        Chromosome contextChromosome = HumanChromosome.fromString(normal.chromosome());
        int minGermlineCoverage =
                contextChromosome.isAllosome() ? config.minGermlineReadContextCoverageAllosome() : config.minGermlineReadContextCoverage();
        if (normal.primaryReadContext().coverage() < minGermlineCoverage) {
            result.add(SoftFilter.MIN_GERMLINE_DEPTH.toString());
        }

        if (Doubles.greaterThan(normalCounter.vaf(), config.maxGermlineVaf())) {
            result.add(SoftFilter.MAX_GERMLINE_VAF.toString());
        }

        double tumorQual = primaryTumor.rawBaseQualityAlt();
        double germlineQual = normal.rawBaseQualityAlt();
        if (Doubles.positive(tumorQual)) {
            if (Doubles.greaterThan(germlineQual / tumorQual, config.maxGermlineRelativeQual())) {
                result.add(SoftFilter.MAX_GERMLINE_REL_RAW_BASE_QUAL.toString());
            }
        }

        // MNV Tests
        if (tier != SageVariantTier.HOTSPOT && normal.isMNV()) {
            if (normal.primaryReadContext().altSupport() != 0) {
                result.add(SoftFilter.MAX_GERMLINE_ALT_SUPPORT.toString());
            }
        }

        return result;
    }

    private boolean skipMinTumorQualTest(@NotNull final SageVariantTier tier, @NotNull final AltContext primaryTumor) {
        return tier.equals(SageVariantTier.HOTSPOT) && primaryTumor.rawDepthAlt() >= config.hotspotMinRawTumorAltSupportToSkipQualCheck()
                && Doubles.greaterOrEqual(primaryTumor.rawVaf(), config.hotspotMinRawTumorVafToSkipQualCheck());
    }

}
