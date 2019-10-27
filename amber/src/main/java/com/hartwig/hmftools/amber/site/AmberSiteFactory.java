package com.hartwig.hmftools.amber.site;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class AmberSiteFactory {

    public static ListMultimap<Chromosome, AmberSite> sites(@NotNull final String vcfFile) throws IOException {
        final ListMultimap<Chromosome, AmberSite> result = ArrayListMultimap.create();

        try (final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFile, new VCFCodec(), false)) {
            for (VariantContext variant : reader.iterator()) {
                if (variant.isNotFiltered()) {
                    if (HumanChromosome.contains(variant.getContig())) {
                        HumanChromosome chromosome = HumanChromosome.fromString(variant.getContig());
                        result.put(chromosome,
                                ImmutableAmberSite.builder()
                                        .chromosome(variant.getContig())
                                        .position(variant.getStart())
                                        .ref(variant.getReference().getBaseString())
                                        .alt(variant.getAlternateAllele(0).getBaseString())
                                        .build());
                    }
                }
            }
        }

        return result;
    }

}
