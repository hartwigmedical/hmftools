package com.hartwig.hmftools.common.variant.enrich;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class NearIndelPonEnrichment implements SomaticEnrichment {

    private static final int DISTANCE = 5;
    private static final String FILTER = "NEAR_INDEL_PON";

    private NearHotspotFilter nearHotspotFilter;

    public NearIndelPonEnrichment(@NotNull final String vcfFile) throws IOException {

        final Multimap<String, VariantHotspot> indelPons = ArrayListMultimap.create();

        try (final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFile, new VCFCodec(), true)) {
            for (final VariantContext context : reader.iterator()) {
                if (context.getType() == VariantContext.Type.INDEL) {
                    indelPons.put(context.getContig(),
                            ImmutableVariantHotspot.builder()
                                    .chromosome(context.getContig())
                                    .position(context.getStart())
                                    .ref(context.getReference().getBaseString())
                                    .alt(context.getAltAlleleWithHighestAlleleCount().getBaseString())
                                    .build());
                }
            }
        }

        nearHotspotFilter = new NearHotspotFilter(DISTANCE, indelPons);
    }

    @NotNull
    @Override
    public ImmutableSomaticVariantImpl.Builder enrich(@NotNull final ImmutableSomaticVariantImpl.Builder builder,
            @NotNull final VariantContext context) {
        if (context.getType() == VariantContext.Type.INDEL & nearHotspotFilter.test(context)) {
            builder.filter(FILTER);
        }

        return builder;
    }
}
