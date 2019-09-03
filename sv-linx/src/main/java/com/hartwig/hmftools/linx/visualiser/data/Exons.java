package com.hartwig.hmftools.linx.visualiser.data;

import java.io.IOException;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.HmfExonRegion;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.region.Strand;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExonFile;

import org.jetbrains.annotations.NotNull;

public class Exons
{

    private static final Comparator<Exon> RANKED = Comparator.comparingInt(Exon::rank);

    @NotNull
    public static List<Exon> geneExons(@NotNull final List<Gene> gene, @NotNull final List<Exon> exons)
    {
        final Set<String> transcripts = gene.stream().map(Gene::transcript).collect(Collectors.toSet());
        return exons.stream()
                .filter(x -> transcripts.contains(x.transcript()) && !x.type().equals(ExonType.DISRUPTED))
                .collect(Collectors.toList());
    }

    @NotNull
    static List<Exon> sortedUpstreamExons(@NotNull final Fusion fusion, @NotNull final List<Exon> exons)
    {
        return exons.stream()
                .filter(x -> x.gene().equals(fusion.geneUp()) & x.transcript().equals(fusion.transcriptUp()))
                .sorted(RANKED)
                .collect(Collectors.toList());
    }

    @NotNull
    static List<Exon> sortedDownstreamExons(@NotNull final Fusion fusion, @NotNull final List<Exon> exons)
    {
        return exons.stream()
                .filter(x -> x.gene().equals(fusion.geneDown()) & x.transcript().equals(fusion.transcriptDown()))
                .sorted(RANKED)
                .collect(Collectors.toList());
    }

    @NotNull
    public static List<Exon> fromHmfTranscript(@NotNull final String sampleId, int clusterId, @NotNull HmfTranscriptRegion transcript)
    {
        final List<Exon> result = Lists.newArrayList();

        for (int i = 0; i < transcript.exome().size(); i++)
        {
            final HmfExonRegion hmfExon = transcript.exome().get(i);
            int rank = transcript.strand().equals(Strand.FORWARD) ? i + 1 : transcript.exome().size() - i;

            Exon exon = ImmutableExon.builder()
                    .type(ExonType.DRIVER)
                    .sampleId(sampleId)
                    .clusterId(clusterId)
                    .gene(transcript.gene())
                    .chromosome(transcript.chromosome())
                    .rank(rank)
                    .start(hmfExon.start())
                    .end(hmfExon.end())
                    .transcript(transcript.transcriptID())
                    .build();

            result.add(exon);
        }

        return result;
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
                .type(ExonType.fromString(file.AnnotationType))
                .sampleId(file.SampleId)
                .clusterId(file.ClusterId)
                .gene(file.Gene)
                .chromosome(file.Chromosome)
                .rank(file.ExonRank)
                .start(file.ExonStart)
                .end(file.ExonEnd)
                .transcript(file.Transcript)
                .build();

    }
}
