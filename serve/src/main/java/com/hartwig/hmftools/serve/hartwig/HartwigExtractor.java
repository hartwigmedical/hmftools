package com.hartwig.hmftools.serve.hartwig;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.hotspot.HotspotGenerator;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class HartwigExtractor {

    private static final Logger LOGGER = LogManager.getLogger(HartwigExtractor.class);

    @NotNull
    private final HotspotGenerator hotspotGenerator;

    public HartwigExtractor(@NotNull final HotspotGenerator hotspotGenerator) {
        this.hotspotGenerator = hotspotGenerator;
    }

    @NotNull
    public Map<HartwigEntry, List<VariantHotspot>> extractFromHartwigEntries(@NotNull List<? extends HartwigEntry> entries) {
        Map<HartwigEntry, List<VariantHotspot>> hotspotsPerEntry = Maps.newHashMap();
        for (HartwigEntry entry : entries) {
            List<VariantHotspot> hotspots = Lists.newArrayList();
            if (!entry.proteinAnnotation().isEmpty()) {
                if (HotspotGenerator.isResolvableProteinAnnotation(entry.proteinAnnotation())) {
                    hotspots = hotspotGenerator.generateHotspots(entry.gene(), entry.transcript(), entry.proteinAnnotation());
                } else {
                    LOGGER.warn("Cannot resolve Hartwig Cohort protein annotation: '{}:p.{}'", entry.gene(), entry.proteinAnnotation());
                }
            }

            VariantHotspot impliedHotspot = toHotspot(entry);
            if (!hotspots.contains(impliedHotspot)) {
                LOGGER.debug("Adding implied hotspot '{}' since it was not resolved from protein annotation '{}'",
                        impliedHotspot,
                        entry.proteinAnnotation());
                hotspots.add(impliedHotspot);
            }
        }

        return hotspotsPerEntry;
    }

    @NotNull
    private static VariantHotspot toHotspot(@NotNull HartwigEntry entry) {
        return ImmutableVariantHotspotImpl.builder()
                .chromosome(entry.chromosome())
                .position(entry.position())
                .ref(entry.ref())
                .alt(entry.alt())
                .build();
    }
}
