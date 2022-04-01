package com.hartwig.hmftools.common.amber;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public final class AmberSiteFactory
{
    public static final Logger LOGGER = LogManager.getLogger(AmberSiteFactory.class);

    @NotNull
    public static AmberSite tumorSite(@NotNull final TumorBAF baseDepth)
    {
        return ImmutableAmberSite.builder()
                .from(baseDepth)
                .snpCheck(false)
                .ref(baseDepth.ref())
                .alt(baseDepth.alt())
                .build();
    }

    @NotNull
    public static AmberSite asSite(@NotNull final BaseDepth baseDepth)
    {
        return ImmutableAmberSite.builder()
                .from(baseDepth)
                .snpCheck(false)
                .ref(baseDepth.ref().toString())
                .alt(baseDepth.alt().toString())
                .build();
    }

    @NotNull
    public static ListMultimap<Chromosome, AmberSite> sites(@NotNull final String vcfFile) throws IOException
    {
        final ListMultimap<Chromosome, AmberSite> result = ArrayListMultimap.create();

        try(final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFile, new VCFCodec(), false))
        {
            for(VariantContext variant : reader.iterator())
            {
                if(variant.isNotFiltered())
                {
                    if(HumanChromosome.contains(variant.getContig()))
                    {
                        HumanChromosome chromosome = HumanChromosome.fromString(variant.getContig());
                        result.put(chromosome,
                                ImmutableAmberSite.builder()
                                        .chromosome(variant.getContig())
                                        .position(variant.getStart())
                                        .ref(variant.getReference().getBaseString())
                                        .alt(variant.getAlternateAllele(0).getBaseString())
                                        .snpCheck(variant.hasAttribute("SNPCHECK"))
                                        .build());
                    }
                }
            }
        }

        LOGGER.info("loaded {} baf loci", result.size());

        return result;
    }
}
