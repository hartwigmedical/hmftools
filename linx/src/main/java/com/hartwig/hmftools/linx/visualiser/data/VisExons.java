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
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExon;

public class VisExons
{
    private static final Comparator<VisGeneExon> RANKED = Comparator.comparingInt(x -> x.ExonRank);

    public static List<VisGeneExon> fusionExons(final Fusion fusion, final List<VisGeneExon> exons)
    {
        final List<VisGeneExon> result = Lists.newArrayList();

        result.addAll(transcriptExons(fusion.transcriptUp(), exons));
        if(!fusion.transcriptUp().equals(fusion.transcriptDown()))
        {
            result.addAll(transcriptExons(fusion.transcriptDown(), exons));
        }

        return result;
    }

    public static List<VisGeneExon> geneExons(final List<Gene> genes, final List<VisGeneExon> exons)
    {
        final List<VisGeneExon> result = Lists.newArrayList();

        for(Gene gene : genes)
        {
            result.addAll(transcriptExons(gene.transcript(), exons));
        }

        return result;
    }

    private static List<VisGeneExon> transcriptExons(final String transcript, final List<VisGeneExon> exons)
    {
        // The same exon might appear more than once so need to remove duplicates
        final Map<Integer,VisGeneExon> dedup = Maps.newLinkedHashMap();

        for(VisGeneExon exon : exons)
        {
            if(exon.Transcript.equals(transcript))
            {
                if(!dedup.containsKey(exon.ExonRank))
                {
                    dedup.put(exon.ExonRank, exon);
                }
            }
        }

        return Lists.newArrayList(dedup.values());
    }

    static List<VisGeneExon> sortedUpstreamExons(final Fusion fusion, final List<VisGeneExon> exons)
    {
        return exons.stream()
                .filter(x -> x.Gene.equals(fusion.geneUp()) & x.Transcript.equals(fusion.transcriptUp()))
                .sorted(RANKED)
                .collect(Collectors.toList());
    }

    static List<VisGeneExon> sortedDownstreamExons(final Fusion fusion, final List<VisGeneExon> exons)
    {
        return exons.stream()
                .filter(x -> x.Gene.equals(fusion.geneDown()) & x.Transcript.equals(fusion.transcriptDown()))
                .sorted(RANKED)
                .collect(Collectors.toList());
    }

    public static List<VisGeneExon> extractExonList(
            final String sampleId, int clusterId, final GeneData geneData, final TranscriptData transcript)
    {
        final List<VisGeneExon> result = transcript.exons().stream()
                .map(x -> new VisGeneExon(
                        sampleId, clusterId, geneData.GeneName, transcript.TransName, geneData.Chromosome, DRIVER,
                        x.Rank, x.Start, x.End))
                .collect(Collectors.toList());

        return result;
    }

    public static List<VisGeneExon> readExons(final String fileName) throws IOException
    {
        return VisGeneExon.read(fileName);
    }
}
