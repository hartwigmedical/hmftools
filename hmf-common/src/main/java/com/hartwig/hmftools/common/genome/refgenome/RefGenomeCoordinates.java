package com.hartwig.hmftools.common.genome.refgenome;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;

public enum RefGenomeCoordinates
{
    COORDS_37(fromResource("lengths.37.tsv"), fromResource("centromeres.37.tsv")),
    COORDS_38(fromResource("lengths.38.tsv"), fromResource("centromeres.38.tsv"));

    @NotNull
    private final Map<Chromosome, Long> lengths;
    @NotNull
    private final Map<Chromosome, Long> centromeres;

    private static final String FIELD_SEPARATOR = "\t";

    RefGenomeCoordinates(@NotNull final Map<Chromosome, Long> lengths, @NotNull final Map<Chromosome, Long> centromeres)
    {
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
    private static Map<Chromosome, Long> fromResource(@NotNull final String resource)
    {
        final InputStream inputStream = RefGenomeCoordinates.class.getResourceAsStream("/refgenome/" + resource);
        return fromLines(new BufferedReader(new InputStreamReader(inputStream)).lines().collect(Collectors.toList()));
    }

    @NotNull
    private static Map<Chromosome, Long> fromLines(@NotNull final List<String> lines)
    {
        final Map<Chromosome, Long> result = Maps.newHashMap();
        for (final String line : lines)
        {
            final String[] values = line.split(FIELD_SEPARATOR);
            result.put(HumanChromosome.fromString(values[0]), Long.valueOf(values[1]));
        }

        return result;
    }
}
