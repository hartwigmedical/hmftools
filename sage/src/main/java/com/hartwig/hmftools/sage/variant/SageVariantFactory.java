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
    public SageVariant create(@NotNull final ReadContextCounter normal) {
        final SageVariantTier tier = tierSelector.tier(normal);
        final Set<String> filters = germlineOnlyFilters(normal);

        return new SageVariant(tier, filters, Collections.singletonList(normal), Collections.emptyList());
    }

    @NotNull
    public SageVariant create(@NotNull final List<ReadContextCounter> normal, @NotNull final List<ReadContextCounter> tumorAltContexts) {

        final ReadContextCounter primaryNormal = normal.get(0);

        final SageVariantTier tier = tierSelector.tier(primaryNormal);
        final SoftFilterConfig softConfig = config.softConfig(tier);

        boolean passingTumor = false;
        final Set<String> allFilters = Sets.newHashSet();
        for (ReadContextCounter tumorAltContext : tumorAltContexts) {
            final Set<String> tumorFilters = pairedFilters(tier, softConfig, primaryNormal, tumorAltContext);
            if (tumorFilters.isEmpty()) {
                passingTumor = true;
            }
            allFilters.addAll(tumorFilters);
        }

        return new SageVariant(tier, passingTumor ? Sets.newHashSet() : allFilters, normal, tumorAltContexts);
    }

    @NotNull
    private Set<String> germlineOnlyFilters(@NotNull final ReadContextCounter germline) {
        final Set<String> result = Sets.newHashSet();

        if (Doubles.lessThan(germline.vaf(), config.minGermlineVaf())) {
            result.add(SageVCF.MIN_GERMLINE_VAF);
        }

        return result;
    }

    @NotNull
    private Set<String> pairedFilters(@NotNull final SageVariantTier tier, @NotNull final SoftFilterConfig config,
            @NotNull final ReadContextCounter normal, @NotNull final ReadContextCounter primaryTumor) {
        Set<String> result = Sets.newHashSet();

        // TUMOR Tests
        final boolean skipTumorTests = skipMinTumorQualTest(tier, primaryTumor);
        if (!skipTumorTests && primaryTumor.tumorQuality() < config.minTumorQual()) {
            result.add(SoftFilter.MIN_TUMOR_QUAL.toString());
        }

        if (!skipTumorTests && Doubles.lessThan(primaryTumor.vaf(), config.minTumorVaf())) {
            result.add(SoftFilter.MIN_TUMOR_VAF.toString());
        }

        // GERMLINE Tests
        Chromosome contextChromosome = HumanChromosome.fromString(normal.chromosome());
        int minGermlineCoverage =
                contextChromosome.isAllosome() ? config.minGermlineReadContextCoverageAllosome() : config.minGermlineReadContextCoverage();
        if (normal.coverage() < minGermlineCoverage) {
            result.add(SoftFilter.MIN_GERMLINE_DEPTH.toString());
        }

        if (Doubles.greaterThan(normal.vaf(), config.maxGermlineVaf())) {
            result.add(SoftFilter.MAX_GERMLINE_VAF.toString());
        }

        double tumorQual = primaryTumor.rawAltBaseQuality();
        double germlineQual = normal.rawAltBaseQuality();
        if (Doubles.positive(tumorQual)) {
            if (Doubles.greaterThan(germlineQual / tumorQual, config.maxGermlineRelativeQual())) {
                result.add(SoftFilter.MAX_GERMLINE_REL_RAW_BASE_QUAL.toString());
            }
        }

        // MNV Tests
        if (tier != SageVariantTier.HOTSPOT && normal.variant().isMNV()) {
            if (normal.altSupport() != 0) {
                result.add(SoftFilter.MAX_GERMLINE_ALT_SUPPORT.toString());
            }
        }

        return result;
    }

    private boolean skipMinTumorQualTest(@NotNull final SageVariantTier tier, @NotNull final ReadContextCounter primaryTumor) {
        return tier.equals(SageVariantTier.HOTSPOT) && primaryTumor.rawAltSupport() >= config.hotspotMinRawTumorAltSupportToSkipQualCheck()
                && Doubles.greaterOrEqual(primaryTumor.rawVaf(), config.hotspotMinRawTumorVafToSkipQualCheck());
    }

}
