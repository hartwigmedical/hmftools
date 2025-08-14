package com.hartwig.hmftools.common.variant;

import java.io.IOException;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public final class VariantHotspotFile
{
    public static ListMultimap<Chromosome,VariantHotspot> readFromVCF(final String fileName) throws IOException
    {
        ListMultimap<Chromosome, VariantHotspot> result = ArrayListMultimap.create();

        VcfFileReader vcfFileReader = new VcfFileReader(fileName);

        for(VariantContext variantContext : vcfFileReader.iterator())
        {
            if(HumanChromosome.contains(variantContext.getContig()))
            {
                result.put(HumanChromosome.fromString(variantContext.getContig()), fromVariantContext(variantContext));
            }
        }

        vcfFileReader.close();

        return result;
    }

    private static VariantHotspot fromVariantContext(@NotNull final VariantContext context)
    {
        return ImmutableVariantHotspotImpl.builder()
                .chromosome(context.getContig())
                .position(context.getStart())
                .ref(context.getReference().getBaseString())
                .alt(context.getAlternateAllele(0).getBaseString())
                .build();
    }
}
