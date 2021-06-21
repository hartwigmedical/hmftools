package com.hartwig.hmftools.linx.visualiser.circos;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegionBuilder;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.linx.visualiser.data.CopyNumberAlteration;
import com.hartwig.hmftools.linx.visualiser.data.Exon;
import com.hartwig.hmftools.linx.visualiser.data.FusedExon;
import com.hartwig.hmftools.linx.visualiser.data.Gene;
import com.hartwig.hmftools.linx.visualiser.data.ImmutableCopyNumberAlteration;
import com.hartwig.hmftools.linx.visualiser.data.ImmutableExon;
import com.hartwig.hmftools.linx.visualiser.data.ImmutableFusedExon;
import com.hartwig.hmftools.linx.visualiser.data.ImmutableGene;
import com.hartwig.hmftools.linx.visualiser.data.ImmutableProteinDomain;
import com.hartwig.hmftools.linx.visualiser.data.ImmutableSegment;
import com.hartwig.hmftools.linx.visualiser.data.ImmutableVisSvData;
import com.hartwig.hmftools.linx.visualiser.data.VisSvData;
import com.hartwig.hmftools.linx.visualiser.data.ProteinDomain;
import com.hartwig.hmftools.linx.visualiser.data.Segment;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class ScalePosition
{
    private static final double MIN_CONTIG_PERCENTAGE = 0.01;
    private static final double GENE_NAME_DISTANCE = 1.5 / 360d;

    private static final Logger LOGGER = LogManager.getLogger(ScalePosition.class);

    private final int totalLength;
    private final Set<GenomePosition> contigLength = Sets.newHashSet();
    private final Map<String, ScaleContig> contigMap = Maps.newHashMap();

    ScalePosition(@NotNull final List<? extends GenomePosition> positions)
    {
        int totalLength = 0;
        final Set<String> contigs = positions.stream().map(GenomePosition::chromosome).collect(Collectors.toSet());
        for (final String contig : contigs)
        {
            final List<Long> contigPositions = positions.stream()
                    .filter(x -> x.chromosome().equals(contig))
                    .map(GenomePosition::position)
                    .collect(Collectors.toList());

            ScaleContig scaleContig = new ScaleContig(contig, contigPositions);
            contigMap.put(contig, scaleContig);
            totalLength += scaleContig.length();
        }

        long minContigDistance = Math.round(MIN_CONTIG_PERCENTAGE * totalLength);

        totalLength = 0;
        for (final ScaleContig scaleContig : contigMap.values())
        {
            if (scaleContig.length() < minContigDistance)
            {
                scaleContig.expand(1d * minContigDistance / scaleContig.length());
            }
            totalLength += scaleContig.length();
            contigLength.add(GenomePositions.create(scaleContig.contig(), scaleContig.length()));
        }

        this.totalLength = totalLength;
    }

    @NotNull
    public Set<GenomePosition> contigLengths()
    {
        return contigLength;
    }

    @NotNull
    public List<Segment> scaleSegments(@NotNull final List<Segment> segments)
    {
        return scale(segments, x -> ImmutableSegment.builder().from(x));
    }


    @NotNull
    public List<CopyNumberAlteration> interpolateAlterations(@NotNull final List<CopyNumberAlteration> segments)
    {
        return segments.stream().map(x -> interpolate(x, y -> ImmutableCopyNumberAlteration.builder().from(x))).collect(Collectors.toList());
    }


    @NotNull
    public List<ProteinDomain> interpolateProteinDomains(@NotNull final List<ProteinDomain> exons)
    {
        return exons.stream().map(x -> interpolate(x, y -> ImmutableProteinDomain.builder().from(x))).collect(Collectors.toList());
    }

    @NotNull
    public List<Exon> interpolateExons(@NotNull final List<Exon> exons)
    {
        return exons.stream().map(x -> interpolate(x, y -> ImmutableExon.builder().from(x))).collect(Collectors.toList());
    }

    @NotNull
    public List<Gene> interpolateGene(@NotNull final List<Gene> genes)
    {
        double geneNameDistance = GENE_NAME_DISTANCE * totalLength;

        return genes.stream().map(x ->
        {
            ScaleContig positionMap = contigMap.get(x.chromosome());
            final int scaledGeneStart = positionMap.interpolate(x.start());
            final int scaledGeneEnd = positionMap.interpolate(x.end());
            int geneNamePosition = (int) Math.round(x.strand() > 0
                    ? scaledGeneStart - geneNameDistance
                    : scaledGeneEnd + geneNameDistance);

            // If gene name won't be written at start of gene because it outside range, move to the end of gene
            if (geneNamePosition < 1)
            {
                geneNamePosition = (int) Math.round(scaledGeneEnd + geneNameDistance);
            }

            if (geneNamePosition > positionMap.length())
            {
                geneNamePosition = (int) Math.round(scaledGeneStart - geneNameDistance);
            }

            return ImmutableGene.builder().from(x).start(scaledGeneStart).end(scaledGeneEnd).namePosition(geneNamePosition).build();
        }).collect(Collectors.toList());
    }


    @NotNull
    public List<FusedExon> scaleFusedExon(@NotNull final List<FusedExon> exons)
    {
        return exons.stream().map(x ->
        {
            ScaleContig positionMap = contigMap.get(x.fusion());
            return scale(x, y -> ImmutableFusedExon.builder()
                    .from(y)
                    .geneStart(positionMap.scale(y.geneStart()))
                    .geneEnd(positionMap.scale(y.geneEnd())), positionMap);
        }).collect(Collectors.toList());
    }

    @NotNull
    private <T extends GenomeRegion> T interpolate(@NotNull final T exon, Function<T, GenomeRegionBuilder<T>> builderFunction)
    {
        final ScaleContig positionMap = contigMap.get(exon.chromosome());
        assert (positionMap != null && !positionMap.isEmpty());

        return builderFunction.apply(exon)
                .start(positionMap.interpolate(exon.start()))
                .end(positionMap.interpolate(exon.end()))
                .build();
    }

    @NotNull
    public List<GenomeRegion> scaleRegions(@NotNull final List<GenomeRegion> regions)
    {
        return regions.stream().map(x -> scale(x, contigMap.get(x.chromosome()))).collect(Collectors.toList());
    }

    @NotNull
    public List<GenomeRegion> interpolateRegions(@NotNull final List<GenomeRegion> regions)
    {
        return regions.stream().map(x -> interpolate(x, contigMap.get(x.chromosome()))).collect(Collectors.toList());
    }

    @NotNull
    public List<VisSvData> scaleLinks(@NotNull final List<VisSvData> links)
    {
        final List<VisSvData> results = Lists.newArrayList();
        for (final VisSvData link : links)
        {
            try
            {
                final ImmutableVisSvData.Builder builder = ImmutableVisSvData.builder().from(link);
                if (link.isValidStart())
                {
                    builder.startPosition(contigMap.get(link.startChromosome()).scale(link.startPosition()));
                }

                if (link.isValidEnd())
                {
                    builder.endPosition(contigMap.get(link.endChromosome()).scale(link.endPosition()));
                }

                results.add(builder.build());
            } catch (Exception e)
            {
                LOGGER.error("Unable to scale link {}", link);
                throw e;
            }
        }

        return results;
    }

    @NotNull
    private static GenomeRegion scale(@NotNull final GenomeRegion region, @NotNull final ScaleContig positionMap)
    {
        return GenomeRegions.create(region.chromosome(), positionMap.scale(region.start()), positionMap.scale(region.end()));
    }

    @NotNull
    private static GenomeRegion interpolate(@NotNull final GenomeRegion region, @NotNull final ScaleContig positionMap)
    {
        return GenomeRegions.create(region.chromosome(), positionMap.interpolate(region.start()), positionMap.interpolate(region.end()));
    }

    @NotNull
    private <T extends GenomeRegion> List<T> scale(@NotNull final List<T> inputs, Function<T, GenomeRegionBuilder<T>> builderFunction)
    {
        return inputs.stream().map(x -> scale(x, builderFunction, contigMap.get(x.chromosome()))).collect(Collectors.toList());
    }

    @NotNull
    private static <T extends GenomeRegion> T scale(@NotNull final T victim, Function<T, GenomeRegionBuilder<T>> builderFunction,
            @NotNull final ScaleContig positionMap)
    {
        return builderFunction.apply(victim)
                .start(positionMap.scale(victim.start()))
                .end(positionMap.scale(victim.end())).build();
    }

}
