package com.hartwig.hmftools.common.refgenome;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Supplier;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;

import org.jetbrains.annotations.NotNull;

public enum RefGenome {
    HG19(fromResource("hg19_len.tsv"), fromResource("hg19_centromere.tsv")),
    HG38(fromResource("hg38_len.tsv"), fromResource("hg38_centromere.tsv"));

    @NotNull
    private final Map<Chromosome, Long> lengths;
    @NotNull
    private final Map<Chromosome, Long> centromeres;

    private static final String FIELD_SEPARATOR = "\t";

    RefGenome(@NotNull final Map<Chromosome, Long> lengths, @NotNull final Map<Chromosome, Long> centromeres) {
        this.lengths = lengths;
        this.centromeres = centromeres;
    }

    @NotNull
    public Map<Chromosome, Long> lengths() {
        return lengths;
    }

    @NotNull
    public Map<Chromosome, Long> centromeres() {
        return centromeres;
    }

    @NotNull
    public static Optional<RefGenome> fromLengths(@NotNull final Map<Chromosome, Long> length) {
        if (length.equals(HG19.lengths())) {
            return Optional.of(HG19);
        }

        if (length.equals(HG38.lengths())) {
            return Optional.of(HG38);
        }

        return Optional.empty();
    }

    @NotNull
    private static Map<Chromosome, Long> fromResource(@NotNull final String resource) {
        final InputStream inputStream = RefGenome.class.getResourceAsStream("/refgenome/" + resource);
        return fromLines(new BufferedReader(new InputStreamReader(inputStream)).lines().collect(Collectors.toList()));
    }

    @NotNull
    private static Map<Chromosome, Long> fromLines(@NotNull final List<String> lines) {
        final Map<Chromosome, Long> result = Maps.newHashMap();
        for (final String line : lines) {
            final String[] values = line.split(FIELD_SEPARATOR);
            result.put(HumanChromosome.fromString(values[0]), Long.valueOf(values[1]));
        }

        return result;
    }
}
