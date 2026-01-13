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
    public static ListMultimap<Chromosome,SimpleVariant> loadHotspotVcf(final String fileName) throws IOException
    {
        ListMultimap<Chromosome,SimpleVariant> result = ArrayListMultimap.create();

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

    private static SimpleVariant fromVariantContext(@NotNull final VariantContext context)
    {
        return new SimpleVariant(
                context.getContig(), context.getStart(), context.getReference().getBaseString(),
                context.getAlternateAllele(0).getBaseString());
    }
}
