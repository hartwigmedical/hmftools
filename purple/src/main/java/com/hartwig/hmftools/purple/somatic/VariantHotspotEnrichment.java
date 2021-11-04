package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.enrich.HotspotEnrichment.DISTANCE;
import static com.hartwig.hmftools.common.variant.enrich.HotspotEnrichment.HOTSPOT_FLAG;
import static com.hartwig.hmftools.common.variant.enrich.HotspotEnrichment.NEAR_HOTSPOT_FLAG;

import java.util.function.Consumer;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.variant.enrich.HotspotEnrichment;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichment;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VariantHotspotEnrichment implements VariantContextEnrichment
{
    // VCF FIELDS
    private static final String HOTSPOT_DESCRIPTION = "Site is at a known hotspot location";
    private static final String NEAR_HOTSPOT_DESCRIPTION = "Variant within " + DISTANCE + " bases of hotspot";

    private final Consumer<VariantContext> mConsumer;
    private final HotspotEnrichment mHotspotEnrichment;

    public VariantHotspotEnrichment(final Multimap<Chromosome, VariantHotspot> hotspots, final Consumer<VariantContext> consumer)
    {
        mConsumer = consumer;
        mHotspotEnrichment = new HotspotEnrichment(hotspots);
    }

    @Override
    public void accept(@NotNull final VariantContext context)
    {
        mConsumer.accept(enrich(context));
    }

    public VariantContext enrich(final VariantContext context)
    {
        VariantContextBuilder builder =
                new VariantContextBuilder(context).attribute(HOTSPOT_FLAG, false).attribute(NEAR_HOTSPOT_FLAG, false);

        if(mHotspotEnrichment.isOnHotspot(context))
        {
            return builder.attribute(HOTSPOT_FLAG, true).make();
        }
        else if(mHotspotEnrichment.isNearHotspot(context))
        {
            return builder.attribute(NEAR_HOTSPOT_FLAG, true).make();
        }

        return builder.make();

    }

    @Override
    public void flush()
    {
        // Not used
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template)
    {
        template.addMetaDataLine(new VCFInfoHeaderLine(HOTSPOT_FLAG, 0, VCFHeaderLineType.Flag, HOTSPOT_DESCRIPTION));
        template.addMetaDataLine(new VCFInfoHeaderLine(NEAR_HOTSPOT_FLAG, 0, VCFHeaderLineType.Flag, NEAR_HOTSPOT_DESCRIPTION));
        return template;
    }
}
