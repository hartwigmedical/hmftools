package com.hartwig.hmftools.linx.visualiser.data;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Comparator;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class FusedExons
{
    private static final String DELIMITER = "\t";
    private static final Comparator<Exon> RANKED = Comparator.comparingInt(Exon::rank);

    public static void write(@NotNull final String fileName, @NotNull final List<FusedExon> fusedExons) throws IOException
    {
        Files.write(new File(fileName).toPath(), toLines(fusedExons));
    }

    @NotNull
    public static List<FusedExon> fusedExons(@NotNull final Fusion fusion, @NotNull final List<Exon> exons)
    {
        final List<FusedExon> result = Lists.newArrayList();

        final List<Exon> upStreamExons =
                exons.stream().filter(x -> x.gene().equals(fusion.geneUp())).sorted(RANKED).collect(Collectors.toList());
        final List<Exon> downStreamExons =
                exons.stream().filter(x -> x.gene().equals(fusion.geneDown())).sorted(RANKED).collect(Collectors.toList());

        final String fusionName = fusion.geneUp() + "_" + fusion.geneDown();
        final long upGeneOffset = upStreamExons.get(0).start();
        final long upGeneStart = Math.abs(upStreamExons.get(0).start() - upGeneOffset) + 1;
        final long upGeneEnd = Math.abs(fusion.positionUp() - upGeneOffset) + 1;

        final ImmutableFusedExon.Builder upFusedExonBuilder = ImmutableFusedExon.builder()
                .sampleId(fusion.sampleId())
                .clusterId(fusion.clusterId())
                .fusion(fusionName)
                .chromosome(fusion.chromosomeUp())
                .gene(fusion.geneUp())
                .geneStart(upGeneStart)
                .geneEnd(upGeneEnd);

        for (final Exon exon : upStreamExons)
        {
            final long exonStart = Math.abs(exon.start() - upGeneOffset) + 1;
            final long exonEnd = Math.abs(exon.end() - upGeneOffset) + 1;

            if (exonStart <= upGeneEnd)
            {
                final FusedExon fusedExon = upFusedExonBuilder
                        .start(exonStart)
                        .end(Math.min(exonEnd, upGeneEnd))
                        .rank(exon.rank())
                        .truncated(exonEnd > upGeneEnd)
                        .build();
                result.add(fusedExon);
            }
        }

        final long downGeneOffset = fusion.positionDown() - upGeneEnd;
        final long downGeneStart = Math.abs(fusion.positionDown() - downGeneOffset) + 1;
        final long downGeneEnd = Math.abs(downStreamExons.get(downStreamExons.size() - 1).end() - downGeneOffset) + 1;

        final ImmutableFusedExon.Builder downFusedExonBuilder = ImmutableFusedExon.builder().from(upFusedExonBuilder.build())
                .chromosome(fusion.chromosomeDown())
                .gene(fusion.geneDown())
                .geneStart(downGeneStart)
                .geneEnd(downGeneEnd);

        for (final Exon exon : downStreamExons)
        {
            final long exonStart = Math.abs(exon.start() - downGeneOffset) + 1;
            final long exonEnd = Math.abs(exon.end() - downGeneOffset) + 1;

            if (exonEnd > downGeneStart)
            {
                final FusedExon fusedExon = downFusedExonBuilder
                        .start(Math.max(exonStart, downGeneStart))
                        .end(exonEnd)
                        .rank(exon.rank())
                        .truncated(exonStart < downGeneStart)
                        .build();
                result.add(fusedExon);
            }

        }

        return result;
    }

    @NotNull
    static List<String> toLines(@NotNull final List<FusedExon> exons)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        exons.stream().map(FusedExons::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(DELIMITER).add("sampleId")
                .add("clusterId")
                .add("fusion")
                .add("gene")
                .add("geneStart")
                .add("geneEnd")
                .add("chromosome")
                .add("start")
                .add("end")
                .add("rank")
                .add("truncated")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final FusedExon exon)
    {
        return new StringJoiner(DELIMITER)
                .add(exon.sampleId())
                .add(String.valueOf(exon.clusterId()))
                .add(String.valueOf(exon.fusion()))
                .add(String.valueOf(exon.gene()))
                .add(String.valueOf(exon.geneStart()))
                .add(String.valueOf(exon.geneEnd()))
                .add(String.valueOf(exon.chromosome()))
                .add(String.valueOf(exon.start()))
                .add(String.valueOf(exon.end()))
                .add(String.valueOf(exon.rank()))
                .add(String.valueOf(exon.truncated()))
                .toString();
    }

}
