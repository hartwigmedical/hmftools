package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.variant.Hotspot.HOTSPOT_DESCRIPTION;
import static com.hartwig.hmftools.common.variant.Hotspot.HOTSPOT_FLAG;
import static com.hartwig.hmftools.common.variant.Hotspot.NEAR_HOTSPOT_FLAG;

import java.util.Collection;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class HotspotEnrichment
{
    // VCF FIELDS
    public static final int HOTSPOT_DISTANCE = 5;
    private static final String NEAR_HOTSPOT_DESCRIPTION = "Variant within " + HOTSPOT_DISTANCE + " bases of hotspot";

    private final Multimap<Chromosome,VariantHotspot> mHotspots;
    private final boolean mEnabled;

    public HotspotEnrichment(final Multimap<Chromosome, VariantHotspot> hotspots, boolean enabled)
    {
        mEnabled = enabled;
        mHotspots = hotspots;
    }

    public void processVariant(final VariantContext variant)
    {
        if(!mEnabled)
            return;

        if(!HumanChromosome.contains(variant.getContig()))
            return;

        if(variant.hasAttribute(HOTSPOT_FLAG) || variant.hasAttribute(NEAR_HOTSPOT_FLAG))
            return;

        final Chromosome chromosome = HumanChromosome.fromString(variant.getContig());
        Collection<VariantHotspot> hotspots = mHotspots.get(chromosome);

        if(hotspots == null)
            return;

        int variantStart = variant.getStart();
        int variantEnd = variantStart + variant.getReference().length() - 1 + HOTSPOT_DISTANCE;

        boolean nearHotspot = false;

        for(VariantHotspot hotspot : hotspots)
        {
            if(variantStart == hotspot.position() && hotspot.ref().equals(variant.getReference().getBaseString())
            && variant.getAlternateAlleles().stream().map(Allele::getBaseString).collect(Collectors.toList()).contains(hotspot.alt()))
            {
                variant.getCommonInfo().putAttribute(HOTSPOT_FLAG, true);
                return;
            }

            int hotspotEnd = hotspot.position() + hotspot.ref().length() - 1 + HOTSPOT_DISTANCE;

            if(positionsOverlap(variantStart, variantEnd, hotspot.position(), hotspotEnd))
            {
                nearHotspot = true;
                // continue checking for exact matches
            }

            if(hotspot.position() > variantEnd)
                break;
        }

        if(nearHotspot)
            variant.getCommonInfo().putAttribute(NEAR_HOTSPOT_FLAG, true);
    }

    public static void enrichHeader(final VCFHeader header)
    {
        header.addMetaDataLine(new VCFInfoHeaderLine(HOTSPOT_FLAG, 0, VCFHeaderLineType.Flag, HOTSPOT_DESCRIPTION));
        header.addMetaDataLine(new VCFInfoHeaderLine(NEAR_HOTSPOT_FLAG, 0, VCFHeaderLineType.Flag, NEAR_HOTSPOT_DESCRIPTION));
    }
}
