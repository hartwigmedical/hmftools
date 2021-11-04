package com.hartwig.hmftools.common.variant.enrich;

import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public class HotspotEnrichment
{
    public static final int DISTANCE = 5;
    public static final String HOTSPOT_FLAG = "HOTSPOT";
    public static final String NEAR_HOTSPOT_FLAG = "NEAR_HOTSPOT";

    private final Multimap<Chromosome, VariantHotspot> mHotspots;

    public HotspotEnrichment(@NotNull final Multimap<Chromosome, VariantHotspot> hotspots)
    {
        this.mHotspots = hotspots;
    }

    @NotNull
    public static Hotspot fromVariant(@NotNull final VariantContext context)
    {
        if(context.getAttributeAsBoolean(HOTSPOT_FLAG, false))
        {
            return Hotspot.HOTSPOT;
        }

        if(context.getAttributeAsBoolean(NEAR_HOTSPOT_FLAG, false))
        {
            return Hotspot.NEAR_HOTSPOT;
        }

        return Hotspot.NON_HOTSPOT;
    }

    public boolean isOnHotspot(@NotNull final VariantContext context)
    {
        if(HumanChromosome.contains(context.getContig()))
        {
            final Chromosome chromosome = HumanChromosome.fromString(context.getContig());
            if(mHotspots.containsKey(chromosome))
            {
                return mHotspots.get(chromosome).stream().anyMatch(x -> exactMatch(x, context));
            }
        }

        return false;
    }

    public boolean isNearHotspot(@NotNull final VariantContext context)
    {
        if(HumanChromosome.contains(context.getContig()))
        {
            final Chromosome chromosome = HumanChromosome.fromString(context.getContig());
            if(mHotspots.containsKey(chromosome))
            {
                return mHotspots.get(chromosome).stream().anyMatch(x -> overlaps(x, context));
            }
        }

        return false;
    }

    private static boolean overlaps(@NotNull final VariantHotspot hotspot, @NotNull final VariantContext variant)
    {
        int variantStart = variant.getStart();
        int variantEnd = variant.getStart() + variant.getReference().length() - 1 + DISTANCE;

        long ponStart = hotspot.position();
        long ponEnd = hotspot.position() + hotspot.ref().length() - 1 + DISTANCE;

        return variantStart <= ponEnd && variantEnd >= ponStart;
    }

    private static boolean exactMatch(@NotNull final VariantHotspot hotspot, @NotNull final VariantContext variant)
    {
        return hotspot.position() == variant.getStart() && hotspot.ref().equals(variant.getReference().getBaseString())
                && variant.getAlternateAlleles().stream().map(Allele::getBaseString).collect(Collectors.toList()).contains(hotspot.alt());
    }
}
