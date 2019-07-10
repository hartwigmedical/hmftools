package com.hartwig.hmftools.linx.visualiser.data;

import static java.util.stream.Collectors.toList;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public class Genes
{

    @NotNull
    public static List<Gene> uniqueGenes(@NotNull final List<Exon> exons) {
        return unique(genes(exons));
    }

    @NotNull
    public static List<Gene> unique(@NotNull final List<Gene> genes)
    {
        final Map<String, Gene> result = Maps.newHashMap();
        for (Gene gene : genes)
        {
            if (result.containsKey(gene.name()))
            {
                final Gene prior = result.get(gene.name());
                if (gene.bases() < prior.bases())
                {
                    continue;
                }
            }
            result.put(gene.name(), gene);
        }

        return Lists.newArrayList(result.values());
    }

    @NotNull
    public static List<Gene> genes(@NotNull final List<Exon> exons)
    {
        final List<Gene> result = Lists.newArrayList();

        final Set<String> transcripts = exons.stream().map(Exon::transcript).collect(Collectors.toSet());
        for (final String transcript : transcripts)
        {
            final List<Exon> transcriptExons = exons.stream().filter(x -> x.transcript().equals(transcript)).sorted().collect(toList());
            final Exon first = transcriptExons.get(0);
            final Exon last = transcriptExons.get(transcriptExons.size() - 1);

            long namePosition = first.rank() <= last.rank() ? first.start() - 1 : last.end() + 1;

            final Gene gene = ImmutableGene.builder()
                    .from(first)
                    .transcript(transcript)
                    .name(first.gene())
                    .end(last.end())
                    .namePosition(namePosition)
                    .build();
            result.add(gene);

        }

        return result;
    }

}
