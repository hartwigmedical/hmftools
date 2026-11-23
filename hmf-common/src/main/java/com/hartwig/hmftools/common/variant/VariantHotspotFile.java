package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.VCF_ZIP_EXTENSION;
import static com.hartwig.hmftools.common.variant.SimpleVariant.loadSimpleVariants;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public final class VariantHotspotFile
{
    public static final String SOMATIC_HOTSPOT_FILE = "somatic_hotspots";
    public static final String SOMATIC_HOTSPOT_FILE_DESC = "Somatic hotspot VCF";

    public static final String GERMLINE_HOTSPOT_FILE = "germline_hotspots";
    public static final String GERMLINE_HOTSPOT_FILE_DESC = "Germline hotspot VCF";

    public static ListMultimap<Chromosome,SimpleVariant> loadHotspotFile(final String fileName) throws IOException
    {
        if(fileName.endsWith(VCF_ZIP_EXTENSION))
            return loadHotspotVcf(fileName);

        ListMultimap<Chromosome,SimpleVariant> result = ArrayListMultimap.create();

        List<SimpleVariant> variants = loadSimpleVariants(fileName);

        for(SimpleVariant variant : variants)
        {
            result.put(HumanChromosome.fromString(variant.chromosome()), variant);
        }

        return result;
    }

    // remove once common resources have switched to only have the TSV file
    private static ListMultimap<Chromosome,SimpleVariant> loadHotspotVcf(final String fileName) throws IOException
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
