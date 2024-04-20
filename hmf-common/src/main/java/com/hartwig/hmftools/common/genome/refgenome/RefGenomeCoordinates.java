package com.hartwig.hmftools.common.genome.refgenome;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

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

    public final Map<Chromosome,Integer> Lengths;
    public final Map<Chromosome,Integer> Centromeres;

    public static final double GENOME_LENGTH_V37 = 3.088e10;
    public static final double GENOME_LENGTH_V38 = 3.096e10;

    RefGenomeCoordinates(@NotNull final Map<Chromosome,Integer> lengths, @NotNull final Map<Chromosome,Integer> centromeres)
    {
        Lengths = lengths;
        Centromeres = centromeres;
    }

    @NotNull
    public Map<Chromosome,Integer> lengths() {
        return Lengths;
    }

    @NotNull
    public Map<Chromosome,Integer> centromeres() {
        return Centromeres;
    }

    public int length(final String chromosome)
    {
        if(!HumanChromosome.contains(chromosome))
            return 0;

        return Lengths.get(HumanChromosome.fromString(chromosome)).intValue();
    }

    public int centromere(final String chromosome)
    {
        if(!HumanChromosome.contains(chromosome))
            return 0;

        return Centromeres.get(HumanChromosome.fromString(chromosome)).intValue();
    }

    private static Map<Chromosome,Integer> fromResource(final String resource)
    {
        final InputStream inputStream = RefGenomeCoordinates.class.getResourceAsStream("/refgenome/" + resource);
        return fromLines(new BufferedReader(new InputStreamReader(inputStream)).lines().collect(Collectors.toList()));
    }

    private static Map<Chromosome,Integer> fromLines(final List<String> lines)
    {
        final Map<Chromosome,Integer> result = Maps.newHashMap();
        for (final String line : lines)
        {
            final String[] values = line.split(TSV_DELIM);
            result.put(HumanChromosome.fromString(values[0]), Integer.valueOf(values[1]));
        }

        return result;
    }
}
