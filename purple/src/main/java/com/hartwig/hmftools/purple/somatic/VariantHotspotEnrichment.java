package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.Hotspot.HOTSPOT_FLAG;
import static com.hartwig.hmftools.common.variant.Hotspot.NEAR_HOTSPOT_FLAG;

import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichment;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VariantHotspotEnrichment
{
    // VCF FIELDS
    public static final int HOTSPOT_DISTANCE = 5;
    private static final String HOTSPOT_DESCRIPTION = "Site is at a known hotspot location";
    private static final String NEAR_HOTSPOT_DESCRIPTION = "Variant within " + HOTSPOT_DISTANCE + " bases of hotspot";

    private final Multimap<Chromosome, VariantHotspot> mHotspots;
    private final boolean mEnabled;

    public VariantHotspotEnrichment(final Multimap<Chromosome, VariantHotspot> hotspots, boolean enabled)
    {
        mEnabled = enabled;
        mHotspots = hotspots;
    }

    public VariantContext processVariant(final VariantContext context)
    {
        if(!mEnabled)
            return context;
        else
            return enrich(context);
    }

    private VariantContext enrich(final VariantContext context)
    {
        VariantContextBuilder builder =
                new VariantContextBuilder(context).attribute(HOTSPOT_FLAG, false).attribute(NEAR_HOTSPOT_FLAG, false);

        if(isOnHotspot(context))
        {
            return builder.attribute(HOTSPOT_FLAG, true).make();
        }
        else if(isNearHotspot(context))
        {
            return builder.attribute(NEAR_HOTSPOT_FLAG, true).make();
        }

        return builder.make();
    }

    public static VCFHeader enrichHeader(@NotNull final VCFHeader template)
    {
        template.addMetaDataLine(new VCFInfoHeaderLine(HOTSPOT_FLAG, 0, VCFHeaderLineType.Flag, HOTSPOT_DESCRIPTION));
        template.addMetaDataLine(new VCFInfoHeaderLine(NEAR_HOTSPOT_FLAG, 0, VCFHeaderLineType.Flag, NEAR_HOTSPOT_DESCRIPTION));
        return template;
    }

    private boolean isOnHotspot(@NotNull final VariantContext context)
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

    private boolean isNearHotspot(@NotNull final VariantContext context)
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
        int variantEnd = variant.getStart() + variant.getReference().length() - 1 + HOTSPOT_DISTANCE;

        long ponStart = hotspot.position();
        long ponEnd = hotspot.position() + hotspot.ref().length() - 1 + HOTSPOT_DISTANCE;

        return variantStart <= ponEnd && variantEnd >= ponStart;
    }

    private static boolean exactMatch(@NotNull final VariantHotspot hotspot, @NotNull final VariantContext variant)
    {
        return hotspot.position() == variant.getStart() && hotspot.ref().equals(variant.getReference().getBaseString())
                && variant.getAlternateAlleles().stream().map(Allele::getBaseString).collect(Collectors.toList()).contains(hotspot.alt());
    }

}
