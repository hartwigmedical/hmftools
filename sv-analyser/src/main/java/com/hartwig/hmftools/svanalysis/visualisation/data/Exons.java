package com.hartwig.hmftools.svanalysis.visualisation.data;

import static com.hartwig.hmftools.svanalysis.visualisation.circos.Span.maxPositionPerChromosome;
import static com.hartwig.hmftools.svanalysis.visualisation.circos.Span.minPositionPerChromosome;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.svanalysis.visualiser.VisGeneExonFile;

import org.jetbrains.annotations.NotNull;

public class Exons
{
    @NotNull
    public static Set<String> genesInSegmentsAndLinks(@NotNull final List<Exon> exons, @NotNull final List<GenomePosition> allPositions)
    {

        final Map<String, Long> minPositionPerChromosome = minPositionPerChromosome(allPositions);
        final Map<String, Long> maxPositionPerChromosome = maxPositionPerChromosome(allPositions);

        final Predicate<Exon> inSegments = exon ->
        {
            long min = minPositionPerChromosome.containsKey(exon.chromosome()) ? minPositionPerChromosome.get(exon.chromosome()) : 0;
            long max = maxPositionPerChromosome.containsKey(exon.chromosome()) ? maxPositionPerChromosome.get(exon.chromosome()) : 0;
            return exon.start() <= max && exon.end() >= min;
        };

        return exons.stream().filter(inSegments).map(Exon::gene).collect(Collectors.toSet());
    }

    @NotNull
    public static List<Exon> readExons(@NotNull final String fileName) throws IOException
    {
        return VisGeneExonFile.read(fileName).stream().map(Exons::fromVis).collect(Collectors.toList());
    }

    @NotNull
    private static Exon fromVis(@NotNull final VisGeneExonFile file)
    {
        return ImmutableExon.builder()
                .sampleId(file.SampleId)
                .clusterId(file.ClusterId)
                .gene(file.Gene)
                .chromosome(file.Chromosome)
                .rank(file.ExonRank)
                .start(file.ExonStart)
                .end(file.ExonEnd)
                .build();

    }
}
