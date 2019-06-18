package com.hartwig.hmftools.linx.visualiser.data;

import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionFactory;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExonFile;

import org.jetbrains.annotations.NotNull;

public class Exons
{

    public static Collection<GenomeRegion> geneSpanPerChromosome(@NotNull final List<Exon> exons)
    {
        final Map<String, GenomeRegion> resultMap = Maps.newHashMap();
        for (Exon exon : exons)
        {
            final String contig = exon.chromosome();

            final GenomeRegion currentGene = resultMap.computeIfAbsent(contig, x -> exon);
            final GenomeRegion newGene =
                    GenomeRegionFactory.create(contig, Math.min(currentGene.start(), exon.start()), Math.max(currentGene.end(), exon
                            .end()));
            resultMap.put(contig, newGene);

        }

        return resultMap.values();
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
