package com.hartwig.hmftools.common.variant.enrich;

import java.io.IOException;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.hotspot.VariantHotspotFile;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public class HotspotEnrichment implements SomaticEnrichment {

    public static final int DISTANCE = 5;
    private static final String HOTSPOT_TAG = "HOTSPOT";

    private final Multimap<Chromosome, VariantHotspot> hotspots;

    @NotNull
    public static HotspotEnrichment fromHotspotsFile(@NotNull String hotspotsFile) throws IOException {
        return new HotspotEnrichment(VariantHotspotFile.read(hotspotsFile));
    }

    @VisibleForTesting
    HotspotEnrichment(@NotNull final Multimap<Chromosome, VariantHotspot> hotspots) {
        this.hotspots = hotspots;
    }

    @NotNull
    @Override
    public ImmutableSomaticVariantImpl.Builder enrich(@NotNull final ImmutableSomaticVariantImpl.Builder builder,
            @NotNull final VariantContext context) {
        if (context.hasAttribute(HOTSPOT_TAG)) {
            return builder.hotspot(Hotspot.HOTSPOT);
        }

        if (isOnHotspot(context)) {
            return builder.hotspot(Hotspot.HOTSPOT);
        }

        if (isNearHotspot(context)) {
            return builder.hotspot(Hotspot.NEAR_HOTSPOT);
        }

        return builder.hotspot(Hotspot.NON_HOTSPOT);
    }

    public boolean isOnHotspot(@NotNull final VariantContext context) {
        if (HumanChromosome.contains(context.getContig())) {
            final Chromosome chromosome = HumanChromosome.fromString(context.getContig());
            if (hotspots.containsKey(chromosome)) {
                return hotspots.get(chromosome).stream().anyMatch(x -> exactMatch(x, context));
            }
        }

        return false;
    }

    public boolean isNearHotspot(@NotNull final VariantContext context) {
        if (HumanChromosome.contains(context.getContig())) {
            final Chromosome chromosome = HumanChromosome.fromString(context.getContig());
            if (hotspots.containsKey(chromosome)) {
                return hotspots.get(chromosome).stream().anyMatch(x -> overlaps(x, context));
            }
        }

        return false;
    }

    private static boolean overlaps(@NotNull final VariantHotspot hotspot, @NotNull final VariantContext variant) {
        int variantStart = variant.getStart();
        int variantEnd = variant.getStart() + variant.getReference().length() - 1 + DISTANCE;

        long ponStart = hotspot.position();
        long ponEnd = hotspot.position() + hotspot.ref().length() - 1 + DISTANCE;

        return variantStart <= ponEnd && variantEnd >= ponStart;
    }

    private static boolean exactMatch(@NotNull final VariantHotspot hotspot, @NotNull final VariantContext variant) {
        return hotspot.position() == variant.getStart() && hotspot.ref().equals(variant.getReference().getBaseString())
                && variant.getAlternateAlleles().stream().map(Allele::getBaseString).collect(Collectors.toList()).contains(hotspot.alt());
    }
}
