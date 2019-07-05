package com.hartwig.hmftools.linx.visualiser.data;

import static java.util.stream.Collectors.toList;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class Genes
{

    @NotNull
    public static List<Gene> downstreamGene(@NotNull final List<Fusion> fusions, @NotNull final List<Exon> exons)
    {
        final List<Gene> result = Lists.newArrayList();

        for (Fusion fusion : fusions)
        {
            final List<Exon> fusionExons = Exons.downstreamExons(fusion, exons);
            if (!fusionExons.isEmpty())
            {
                result.add(FusedExons.downGeneRegion(fusion, fusionExons.get(fusionExons.size() - 1)));
            }

        }

        return result;
    }

    @NotNull
    public static List<Gene> genes(@NotNull final List<Exon> exons)
    {
        final List<Gene> result = Lists.newArrayList();

        final Set<String> genes = exons.stream().map(Exon::gene).collect(Collectors.toSet());
        for (final String geneName : genes)
        {

            final List<Exon> geneExons = exons.stream().filter(x -> x.gene().equals(geneName)).sorted().collect(toList());
            final Exon first = geneExons.get(0);
            final Exon last = geneExons.get(geneExons.size() - 1);

            long namePosition = first.rank() <= last.rank() ? first.start() - 1 : last.end() + 1;

            final Gene gene = ImmutableGene.builder().from(first).name(first.gene()).end(last.end()).namePosition(namePosition).build();
            result.add(gene);

        }

        return result;
    }

}
