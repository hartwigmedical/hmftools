package com.hartwig.hmftools.linx.visualiser.data;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.linx.visualiser.data.VisExons.sortedDownstreamExons;
import static com.hartwig.hmftools.linx.visualiser.data.VisExons.sortedUpstreamExons;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.linx.visualiser.file.VisFusion;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExon;

public class FusedExons
{
    public static void write(final String fileName, final List<FusedExon> fusedExons) throws IOException
    {
        Files.write(new File(fileName).toPath(), toLines(fusedExons));
    }

    public static List<FusedExon> fusedExons(final VisFusion fusion, final List<VisGeneExon> exons)
    {
        final List<FusedExon> result = Lists.newArrayList();

        final List<VisGeneExon> fusionExons = VisExons.fusionExons(fusion, exons);
        final List<VisGeneExon> upStreamExons = sortedUpstreamExons(fusion, fusionExons);
        final List<VisGeneExon> downStreamExons = sortedDownstreamExons(fusion, fusionExons);

        if(upStreamExons.isEmpty() || downStreamExons.isEmpty())
            return result;

        final VisGeneExon firstUpExon = upStreamExons.get(0);
        final GenomeRegion upGeneRegion = formGeneRegion(fusion, firstUpExon, true);
        final GenomeRegion convertedUpGeneRegion = convertRegion(fusion.StrandUp, upGeneRegion, upGeneRegion);

        final ImmutableFusedExon.Builder upFusedExonBuilder = ImmutableFusedExon.builder()
                .sampleId(fusion.SampleId)
                .clusterId(fusion.ClusterId)
                .fusion(fusion.name())
                .chromosome(fusion.ChrUp)
                .unadjustedGeneStart(upGeneRegion.start())
                .gene(fusion.GeneNameUp)
                .geneStart(convertedUpGeneRegion.start())
                .geneEnd(convertedUpGeneRegion.end())
                .isUpstream(true)
                .transcript(fusion.TranscriptUp);

        boolean hasUpExons = false;

        for(final VisGeneExon exon : upStreamExons)
        {
            final GenomeRegion convertedExon = convertRegion(fusion.StrandUp, upGeneRegion, exon);

            if(upGeneRegion.overlaps(exon))
            {
                final FusedExon fusedExon = upFusedExonBuilder
                        .start(convertedExon.start())
                        .end(convertedExon.end())
                        .rank(exon.ExonRank)
                        .skipped(exon.ExonRank > fusion.FusedExonUp)
                        .build();
                result.add(fusedExon);
                hasUpExons = true;
            }
        }

        final VisGeneExon finalDownExon = downStreamExons.get(downStreamExons.size() - 1);
        final GenomeRegion downGeneRegion = formGeneRegion(fusion, finalDownExon, false);
        final GenomeRegion convertedDownGeneRegion = convertRegion(fusion.StrandDown, downGeneRegion, downGeneRegion);

        final ImmutableFusedExon.Builder downFusedExonBuilder = ImmutableFusedExon.builder()
                .sampleId(fusion.SampleId)
                .clusterId(fusion.ClusterId)
                .fusion(fusion.name())
                .chromosome(fusion.ChrDown)
                .unadjustedGeneStart(fusion.PosDown)
                .gene(fusion.GeneNameDown)
                .geneStart(convertedDownGeneRegion.start() + convertedUpGeneRegion.end())
                .geneEnd(convertedDownGeneRegion.end() + convertedUpGeneRegion.end())
                .isUpstream(false)
                .transcript(fusion.TranscriptDown);

        boolean intronicToExonicFusion = fusion.RegionTypeUp.equals("Intronic") && fusion.RegionTypeDown.equals("Exonic");
        boolean hasDownExons = false;

        for(int i = 0; i < downStreamExons.size(); i++)
        {
            final VisGeneExon exon = downStreamExons.get(i);
            final GenomeRegion convertedExon = convertRegion(fusion.StrandDown, downGeneRegion, exon);

            if(downGeneRegion.overlaps(exon))
            {
                final FusedExon fusedExon = downFusedExonBuilder
                        .start(convertedExon.start() + convertedUpGeneRegion.end())
                        .end(convertedExon.end() + convertedUpGeneRegion.end())
                        .rank(exon.ExonRank)
                        .skipped(exon.ExonRank < fusion.FusedExonDown || (i == 0 && intronicToExonicFusion))
                        .build();
                result.add(fusedExon);
                hasDownExons = true;
            }
        }

        if(!hasUpExons || !hasDownExons)
            return Collections.emptyList();

        return result;
    }

    private static Gene formGeneRegion(final VisFusion fusion, final VisGeneExon firstExon, boolean isUpstreamGene)
    {
        // note the first exon is by rank (ie 1) not by lowest position
        ImmutableGene.Builder builder = ImmutableGene.builder();

        builder.type(firstExon.AnnotationType);
        builder.namePosition(0);
        builder.chromosome(firstExon.Chromosome);

        if(isUpstreamGene)
        {
            builder.name(fusion.GeneNameUp);
            builder.transcript(fusion.TranscriptUp);
            builder.strand(fusion.StrandUp);

            if(fusion.StrandUp > 0)
            {
                  builder.start(firstExon.ExonStart);
                  builder.end(fusion.PosUp);
            }
            else
            {
                builder.start(fusion.PosUp);
                builder.end(firstExon.ExonEnd);
            }
        }
        else
        {
            builder.name(fusion.GeneNameDown);
            builder.transcript(fusion.TranscriptDown);
            builder.strand(fusion.StrandDown);

            if(fusion.StrandDown > 0)
            {
                builder.start(fusion.PosDown);
                builder.end(firstExon.ExonEnd);
            }
            else
            {
                builder.start(firstExon.ExonStart);
                builder.end(fusion.PosDown);
            }
        }

        return builder.build();
    }

    protected static GenomeRegion convertRegion(int strand, final GenomeRegion reference, final GenomeRegion region)
    {
        final int start;
        final int end;
        if(strand < 0)
        {
            start = reference.end() - region.end();
            end = reference.end() - region.start();
        }
        else
        {
            start = region.start() - reference.start();
            end = Math.min(reference.end(), region.end()) - reference.start();
        }

        return GenomeRegions.create(region.chromosome(), Math.max(0, start), end);
    }

    private static List<String> toLines(final List<FusedExon> exons)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        exons.stream().map(FusedExons::toString).forEach(lines::add);
        return lines;
    }

    private static String header()
    {
        return new StringJoiner(TSV_DELIM).add("sampleId")
                .add("clusterId")
                .add("fusion")
                .add("gene")
                .add("geneStart")
                .add("geneEnd")
                .add("isUpstream")
                .add("chromosome")
                .add("start")
                .add("end")
                .add("rank")
                .add("skipped")
                .add("transcript")
                .toString();
    }

    private static String toString(final FusedExon exon)
    {
        return new StringJoiner(TSV_DELIM)
                .add(exon.sampleId())
                .add(String.valueOf(exon.clusterId()))
                .add(exon.fusion())
                .add(exon.gene())
                .add(String.valueOf(exon.geneStart()))
                .add(String.valueOf(exon.geneEnd()))
                .add(String.valueOf(exon.isUpstream()))
                .add(exon.chromosome())
                .add(String.valueOf(exon.start()))
                .add(String.valueOf(exon.end()))
                .add(String.valueOf(exon.rank()))
                .add(String.valueOf(exon.skipped()))
                .add(exon.transcript())
                .toString();
    }
}
