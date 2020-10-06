package com.hartwig.hmftools.common.variant.hotspot;

import java.io.IOException;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public final class VariantHotspotFile {

    private VariantHotspotFile() {
    }

    @NotNull
    public static ListMultimap<Chromosome, VariantHotspot> readFromVCF(@NotNull final String fileName) throws IOException {
        ListMultimap<Chromosome, VariantHotspot> result = ArrayListMultimap.create();

        try (final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(fileName,
                new VCFCodec(),
                false)) {
            for (VariantContext variantContext : reader.iterator()) {
                if (HumanChromosome.contains(variantContext.getContig())) {
                    result.put(HumanChromosome.fromString(variantContext.getContig()), fromVariantContext(variantContext));
                }
            }
        }

        return result;
    }


    @NotNull
    private static VariantHotspot fromVariantContext(@NotNull final VariantContext context) {
        return ImmutableVariantHotspotImpl.builder()
                .chromosome(context.getContig())
                .position(context.getStart())
                .ref(context.getReference().getBaseString())
                .alt(context.getAlternateAllele(0).getBaseString())
                .build();
    }
}
