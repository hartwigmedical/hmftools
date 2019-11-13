package com.hartwig.hmftools.common.hotspot;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public final class VariantHotspotFile {

    private static final String DELIMITER = "\t";

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
    public static ListMultimap<Chromosome, VariantHotspot> read(@NotNull final String fileName) throws IOException {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    @NotNull
    private static ListMultimap<Chromosome, VariantHotspot> fromLines(@NotNull List<String> lines) {
        ListMultimap<Chromosome, VariantHotspot> result = ArrayListMultimap.create();
        for (String line : lines) {
            VariantHotspot position = fromString(line);
            if (HumanChromosome.contains(position.chromosome())) {
                result.put(HumanChromosome.fromString(position.chromosome()), position);
            }
        }

        return result;
    }

    @NotNull
    private static VariantHotspot fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        return ImmutableVariantHotspotImpl.builder()
                .chromosome(values[0])
                .position(Long.parseLong(values[1]))
                .ref(values[2])
                .alt(values[3])
                .build();
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
