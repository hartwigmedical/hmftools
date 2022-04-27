package com.hartwig.hmftools.linx.visualiser.data;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.linx.visualiser.file.VisGeneAnnotationType.EXON_LOST;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExon;

import org.jetbrains.annotations.NotNull;

public class Genes
{
    public static List<Gene> uniqueGenes(final List<VisGeneExon> exons)
    {
        return unique(genes(exons));
    }

    public static List<Gene> unique(final List<Gene> genes)
    {
        final Map<String, Gene> result = Maps.newHashMap();
        for(Gene gene : genes)
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

    public static List<Gene> genes(final List<VisGeneExon> exons)
    {
        final List<Gene> result = Lists.newArrayList();

        final Set<String> transcripts = exons.stream()
                .filter(x -> x.AnnotationType != EXON_LOST)
                .map(x -> x.Transcript)
                .collect(Collectors.toSet());

        for(final String transcript : transcripts)
        {
            final List<VisGeneExon> transcriptExons = exons.stream().filter(x -> x.Transcript.equals(transcript)).sorted().collect(toList());
            final VisGeneExon first = transcriptExons.get(0);
            final VisGeneExon last = transcriptExons.get(transcriptExons.size() - 1);

            long namePosition = first.ExonRank <= last.ExonRank ? first.ExonStart : last.ExonEnd;

            final Gene gene = ImmutableGene.builder()
                    .from(first)
                    .type(first.AnnotationType)
                    .transcript(first.Transcript)
                    .name(first.Gene)
                    .end(last.ExonEnd)
                    .namePosition(namePosition)
                    .strand(first.ExonRank <= last.ExonRank ? 1 : -1)
                    .build();
            result.add(gene);
        }

        return result;
    }

}
