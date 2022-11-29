package com.hartwig.hmftools.linx.visualiser.circos;

import static com.hartwig.hmftools.linx.visualiser.SvVisualiser.VIS_LOGGER;
import static com.hartwig.hmftools.linx.visualiser.file.VisSegment.from;

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
import com.hartwig.hmftools.linx.visualiser.data.FusedExon;
import com.hartwig.hmftools.linx.visualiser.data.Gene;
import com.hartwig.hmftools.linx.visualiser.data.ImmutableFusedExon;
import com.hartwig.hmftools.linx.visualiser.data.ImmutableGene;
import com.hartwig.hmftools.linx.visualiser.file.VisCopyNumber;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExon;
import com.hartwig.hmftools.linx.visualiser.file.VisProteinDomain;
import com.hartwig.hmftools.linx.visualiser.file.VisSegment;
import com.hartwig.hmftools.linx.visualiser.file.VisSvData;

import org.jetbrains.annotations.NotNull;

class ScalePosition
{
    private static final double MIN_CONTIG_PERCENTAGE = 0.01;
    private static final double GENE_NAME_DISTANCE = 1.5 / 360d;

    private final int mTotalLength;
    private final Set<GenomePosition> mContigLength = Sets.newHashSet();
    private final Map<String, ScaleContig> mContigMap = Maps.newHashMap();

    ScalePosition(final List<? extends GenomePosition> positions)
    {
        int totalLength = 0;
        final Set<String> contigs = positions.stream().map(GenomePosition::chromosome).collect(Collectors.toSet());
        for (final String contig : contigs)
        {
            final List<Integer> contigPositions = positions.stream()
                    .filter(x -> x.chromosome().equals(contig))
                    .map(GenomePosition::position)
                    .collect(Collectors.toList());

            ScaleContig scaleContig = new ScaleContig(contig, contigPositions);
            mContigMap.put(contig, scaleContig);
            totalLength += scaleContig.length();
        }

        long minContigDistance = Math.round(MIN_CONTIG_PERCENTAGE * totalLength);

        totalLength = 0;
        for (final ScaleContig scaleContig : mContigMap.values())
        {
            if(scaleContig.length() < minContigDistance)
            {
                scaleContig.expand(1d * minContigDistance / scaleContig.length());
            }
            totalLength += scaleContig.length();
            mContigLength.add(GenomePositions.create(scaleContig.contig(), scaleContig.length()));
        }

        mTotalLength = totalLength;
    }

    public Set<GenomePosition> contigLengths() { return mContigLength; }

    public List<VisSegment> scaleSegments(final List<VisSegment> segments)
    {
        List<VisSegment> results = Lists.newArrayList();

        for(VisSegment segment : segments)
        {
            VisSegment newSegment = from(segment);
            ScaleContig positionMap = mContigMap.get(segment.Chromosome);
            newSegment.setStart(positionMap.scale(segment.start()));
            newSegment.setEnd(positionMap.scale(segment.end()));
            results.add(newSegment);
        }

        return results;
    }

    public List<VisCopyNumber> interpolateCopyNumbers(final List<VisCopyNumber> segments)
    {
        List<VisCopyNumber> results = Lists.newArrayList();

        for(VisCopyNumber cnData : segments)
        {
            VisCopyNumber newCnData = VisCopyNumber.from(cnData);
            ScaleContig positionMap = mContigMap.get(cnData.Chromosome);
            newCnData.Start = positionMap.scale(cnData.start());
            newCnData.End = positionMap.scale(cnData.end());
            results.add(newCnData);
        }

        return results;
    }

    public List<VisProteinDomain> interpolateProteinDomains(final List<VisProteinDomain> exons)
    {
        List<VisProteinDomain> results = Lists.newArrayList();
        for(VisProteinDomain exon : exons)
        {
            VisProteinDomain newExon = VisProteinDomain.from(exon);

            final ScaleContig positionMap = mContigMap.get(exon.chromosome());
            newExon.Start = positionMap.interpolate(exon.start());
            newExon.End = positionMap.interpolate(exon.end());
            results.add(newExon);
        }
        return results;
    }

    public List<VisGeneExon> interpolateExons(final List<VisGeneExon> exons)
    {
        List<VisGeneExon> interpolatedExons = Lists.newArrayList();

        for(VisGeneExon exon : exons)
        {
            final ScaleContig positionMap = mContigMap.get(exon.chromosome());

            interpolatedExons.add(new VisGeneExon(
                    exon.SampleId, exon.ClusterId, exon.Gene, exon.Transcript, exon.Chromosome, exon.AnnotationType, exon.ExonRank,
                    positionMap.interpolate(exon.ExonStart), positionMap.interpolate(exon.ExonEnd)));
        }

        return interpolatedExons;
    }

    public List<Gene> interpolateGene(final List<Gene> genes)
    {
        double geneNameDistance = GENE_NAME_DISTANCE * mTotalLength;

        return genes.stream().map(x ->
        {
            ScaleContig positionMap = mContigMap.get(x.chromosome());
            final int scaledGeneStart = positionMap.interpolate(x.start());
            final int scaledGeneEnd = positionMap.interpolate(x.end());
            int geneNamePosition = (int) Math.round(x.strand() > 0
                    ? scaledGeneStart - geneNameDistance
                    : scaledGeneEnd + geneNameDistance);

            // If gene name won't be written at start of gene because it outside range, move to the end of gene
            if(geneNamePosition < 1)
            {
                geneNamePosition = (int) Math.round(scaledGeneEnd + geneNameDistance);
            }

            if(geneNamePosition > positionMap.length())
            {
                geneNamePosition = (int) Math.round(scaledGeneStart - geneNameDistance);
            }

            return ImmutableGene.builder().from(x).start(scaledGeneStart).end(scaledGeneEnd).namePosition(geneNamePosition).build();
        }).collect(Collectors.toList());
    }


    @NotNull
    public List<FusedExon> scaleFusedExon(final List<FusedExon> exons)
    {
        return exons.stream().map(x ->
        {
            ScaleContig positionMap = mContigMap.get(x.fusion());
            return scale(x, y -> ImmutableFusedExon.builder()
                    .from(y)
                    .geneStart(positionMap.scale(y.geneStart()))
                    .geneEnd(positionMap.scale(y.geneEnd())), positionMap);
        }).collect(Collectors.toList());
    }

    private <T extends GenomeRegion> T interpolate(final T exon, Function<T, GenomeRegionBuilder<T>> builderFunction)
    {
        final ScaleContig positionMap = mContigMap.get(exon.chromosome());
        assert (positionMap != null && !positionMap.isEmpty());

        return builderFunction.apply(exon)
                .start(positionMap.interpolate(exon.start()))
                .end(positionMap.interpolate(exon.end()))
                .build();
    }

    public List<GenomeRegion> interpolateRegions(final List<GenomeRegion> regions)
    {
        return regions.stream().map(x -> interpolate(x, mContigMap.get(x.chromosome()))).collect(Collectors.toList());
    }

    public List<VisSvData> scaleLinks(final List<VisSvData> links)
    {
        final List<VisSvData> results = Lists.newArrayList();
        for (final VisSvData link : links)
        {
            try
            {
                VisSvData newData = VisSvData.from(link);

                if(link.isValidStart())
                {
                    newData.PosStart = mContigMap.get(link.ChrStart).scale(link.PosStart);
                }

                if(link.isValidEnd())
                {
                    newData.PosEnd = mContigMap.get(link.ChrEnd).scale(link.PosEnd);
                }

                results.add(newData);
            }
            catch (Exception e)
            {
                VIS_LOGGER.error("Unable to scale link {}", link);
                throw e;
            }
        }

        return results;
    }

    @NotNull
    private static GenomeRegion interpolate(final GenomeRegion region, final ScaleContig positionMap)
    {
        return GenomeRegions.create(region.chromosome(), positionMap.interpolate(region.start()), positionMap.interpolate(region.end()));
    }

    @NotNull
    private static <T extends GenomeRegion> T scale(final T victim, Function<T, GenomeRegionBuilder<T>> builderFunction,
            final ScaleContig positionMap)
    {
        return builderFunction.apply(victim)
                .start(positionMap.scale(victim.start()))
                .end(positionMap.scale(victim.end())).build();
    }

}
