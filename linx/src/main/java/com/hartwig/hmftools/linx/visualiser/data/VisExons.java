package com.hartwig.hmftools.linx.visualiser.data;

import static com.hartwig.hmftools.linx.visualiser.file.VisGeneAnnotationType.DRIVER;

import java.io.IOException;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneAnnotationType;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExonFile;

import org.jetbrains.annotations.NotNull;

public class VisExons
{
    private static final Comparator<Exon> RANKED = Comparator.comparingInt(Exon::rank);

    @NotNull
    public static List<Exon> fusionExons(@NotNull final Fusion fusion, @NotNull final List<Exon> exons)
    {
        final List<Exon> result = Lists.newArrayList();

        result.addAll(transcriptExons(fusion.transcriptUp(), exons));
        if (!fusion.transcriptUp().equals(fusion.transcriptDown()))
        {
            result.addAll(transcriptExons(fusion.transcriptDown(), exons));
        }

        return result;
    }

    @NotNull
    public static List<Exon> geneExons(@NotNull final List<Gene> genes, @NotNull final List<Exon> exons)
    {
        final List<Exon> result = Lists.newArrayList();

        for (Gene gene : genes)
        {
            result.addAll(transcriptExons(gene.transcript(), exons));
        }

        return result;
    }

    @NotNull
    private static List<Exon> transcriptExons(@NotNull final String transcript, @NotNull final List<Exon> exons)
    {
        // The same exon might appear more than once so need to remove duplicates
        final Map<Integer, Exon> dedup = Maps.newLinkedHashMap();

        for (Exon exon : exons)
        {
            if (exon.transcript().equals(transcript))
            {
                if (!dedup.containsKey(exon.rank()))
                {
                    dedup.put(exon.rank(), exon);
                }
            }
        }

        return Lists.newArrayList(dedup.values());
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
    public static List<Exon> extractExonList(@NotNull final String sampleId, int clusterId, @NotNull HmfTranscriptRegion transcript)
    {
        final List<Exon> result = Lists.newArrayList();

        for (int i = 0; i < transcript.exome().size(); i++)
        {
            final HmfExonRegion hmfExon = transcript.exome().get(i);
            int rank = transcript.strand().equals(Strand.FORWARD) ? i + 1 : transcript.exome().size() - i;

            Exon exon = ImmutableExon.builder()
                    .type(DRIVER)
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
    public static List<Exon> extractExonList(final String sampleId, int clusterId, final GeneData geneData, final TranscriptData transcript)
    {
        final List<Exon> result = Lists.newArrayList();

        for (int i = 0; i < transcript.exons().size(); i++)
        {
            ExonData exonData = transcript.exons().get(i);

            Exon exon = ImmutableExon.builder()
                    .type(DRIVER)
                    .sampleId(sampleId)
                    .clusterId(clusterId)
                    .gene(geneData.GeneName)
                    .chromosome(geneData.Chromosome)
                    .rank(exonData.Rank)
                    .start(exonData.Start)
                    .end(exonData.End)
                    .transcript(transcript.TransName)
                    .build();

            result.add(exon);
        }

        return result;
    }

    @NotNull
    public static List<Exon> readExons(@NotNull final String fileName) throws IOException
    {
        return VisGeneExonFile.read(fileName).stream().map(VisExons::fromVis).collect(Collectors.toList());
    }

    @NotNull
    private static Exon fromVis(@NotNull final VisGeneExonFile file)
    {
        return ImmutableExon.builder()
                .type(VisGeneAnnotationType.valueOf(file.AnnotationType))
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
