package com.hartwig.hmftools.linx.visualiser.data;

import static com.hartwig.hmftools.linx.visualiser.data.Exons.downstreamExons;
import static com.hartwig.hmftools.linx.visualiser.data.Exons.upstreamExons;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegions;

import org.jetbrains.annotations.NotNull;

public class FusedExons
{
    private static final String DELIMITER = "\t";

    public static void write(@NotNull final String fileName, @NotNull final List<FusedExon> fusedExons) throws IOException
    {
        Files.write(new File(fileName).toPath(), toLines(fusedExons));
    }

    @NotNull
    public static List<FusedExon> fusedExons(@NotNull final Fusion fusion, @NotNull final List<Exon> exons)
    {
        final List<FusedExon> result = Lists.newArrayList();

        final List<Exon> upStreamExons = upstreamExons(fusion, exons);
        final List<Exon> downStreamExons = downstreamExons(fusion, exons);
        if (upStreamExons.isEmpty() || downStreamExons.isEmpty())
        {
            return result;
        }

        final Exon firstUpExon = upStreamExons.get(0);
        final GenomeRegion upGeneRegion = upGeneRegion(fusion, firstUpExon);
        final GenomeRegion convertedUpGeneRegion = convertRegion(fusion, upGeneRegion, upGeneRegion);

        final ImmutableFusedExon.Builder upFusedExonBuilder = ImmutableFusedExon.builder()
                .sampleId(fusion.sampleId())
                .clusterId(fusion.clusterId())
                .fusion(fusion.name())
                .chromosome(fusion.chromosomeUp())
                .unadjustedGeneStart(upGeneRegion.start())
                .gene(fusion.geneUp())
                .geneStart(convertedUpGeneRegion.start())
                .geneEnd(convertedUpGeneRegion.end());

        for (final Exon exon : upStreamExons)
        {
            final GenomeRegion convertedExon = convertRegion(fusion, upGeneRegion, exon);

            if (exon.start() <= upGeneRegion.end())
            {
                final FusedExon fusedExon = upFusedExonBuilder
                        .start(convertedExon.start())
                        .end(convertedExon.end())
                        .rank(exon.rank())
                        .skipped(false) // TODO: Check exon skipped field
                        .build();
                result.add(fusedExon);
            }
        }

        final Exon finalDownExon = downStreamExons.get(downStreamExons.size() - 1);
        final GenomeRegion downGeneRegion = downGeneRegion(fusion, finalDownExon);
        final GenomeRegion convertedDownGeneRegion = convertRegion(fusion, downGeneRegion, downGeneRegion);

        final ImmutableFusedExon.Builder downFusedExonBuilder = ImmutableFusedExon.builder().from(upFusedExonBuilder.build())
                .chromosome(fusion.chromosomeDown())
                .unadjustedGeneStart(fusion.positionDown())
                .gene(fusion.geneDown())
                .geneStart(convertedDownGeneRegion.start() + convertedUpGeneRegion.end())
                .geneEnd(convertedDownGeneRegion.end() + convertedUpGeneRegion.end());

        boolean intronicToExonicFusion = fusion.regionTypeUp().equals("Intronic") && fusion.regionTypeDown().equals("Exonic");

        for (int i = 0; i < downStreamExons.size(); i++)
        {
            final Exon exon = downStreamExons.get(i);
            final GenomeRegion convertedExon = convertRegion(fusion, downGeneRegion, exon);

            if (exon.end() > downGeneRegion.start())
            {
                final FusedExon fusedExon = downFusedExonBuilder
                        .start(convertedExon.start() + convertedUpGeneRegion.end())
                        .end(convertedExon.end() + convertedUpGeneRegion.end())
                        .rank(exon.rank())
                        .skipped(exon.rank() == 1 || (i == 0 && intronicToExonicFusion)) // TODO: Check exon skipped field
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
                .add("skipped")
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
                .add(String.valueOf(exon.skipped()))
                .toString();
    }

    @NotNull
    private static GenomeRegion upGeneRegion(@NotNull final Fusion fusion, @NotNull final Exon firstUpGene)
    {
        return fusion.strandUp() < 0
                ? GenomeRegions.create(firstUpGene.chromosome(), fusion.positionUp(), firstUpGene.end())
                : GenomeRegions.create(firstUpGene.chromosome(), firstUpGene.start(), fusion.positionUp());
    }

    @NotNull
    static Gene downGeneRegion(@NotNull final Fusion fusion, @NotNull final Exon finalDownGene)
    {
        return fusion.strandUp() < 0 ?
                ImmutableGene.builder()
                        .chromosome(finalDownGene.chromosome())
                        .start(finalDownGene.start())
                        .end(fusion.positionDown())
                        .namePosition(fusion.positionDown() + 1)
                        .name(fusion.geneDown())
                        .build() :
                ImmutableGene.builder()
                        .chromosome(finalDownGene.chromosome())
                        .start(fusion.positionDown())
                        .end(finalDownGene.end())
                        .namePosition(fusion.positionDown() - 1)
                        .name(fusion.geneDown())
                        .build();
    }

    @NotNull
    static GenomeRegion convertRegion(@NotNull final Fusion fusion, @NotNull final GenomeRegion gene, @NotNull final GenomeRegion region)
    {

        final long start;
        final long end;
        if (fusion.strandUp() < 0)
        {
            start = gene.end() - region.end();
            end = gene.end() - region.start();
        }
        else
        {
            start = region.start() - gene.start();
            end = region.end() - gene.start();
        }

        return GenomeRegions.create(region.chromosome(), start, end);
    }
}
