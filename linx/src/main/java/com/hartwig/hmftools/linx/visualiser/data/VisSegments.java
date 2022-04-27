package com.hartwig.hmftools.linx.visualiser.data;

import static com.hartwig.hmftools.linx.visualiser.circos.Span.maxPositionPerChromosome;
import static com.hartwig.hmftools.linx.visualiser.circos.Span.minPositionPerChromosome;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.linx.visualiser.circos.SegmentTerminal;
import com.hartwig.hmftools.linx.visualiser.circos.Span;
import com.hartwig.hmftools.linx.visualiser.file.VisSegment;
import com.hartwig.hmftools.linx.visualiser.file.VisSvData;

import org.jetbrains.annotations.NotNull;

public class VisSegments
{
    public static double segmentPloidyBefore(final int track, final GenomePosition position, final List<VisSegment> segments)
    {
        return segments.stream()
                .filter(x -> x.Track < track)
                .filter(x -> x.chromosome().equals(position.chromosome()) && (x.start() == position.position()
                        || x.end() == position.position()))
                .mapToDouble(x -> x.LinkPloidy)
                .sum();
    }

    @NotNull
    public static VisSegment entireChromosome(final String sampleId, final String chromosome, final RefGenomeCoordinates refGenCoords)
    {
        VisSegment segment = new VisSegment(
                sampleId, -1, -1, chromosome, "1", String.valueOf(refGenCoords.length(chromosome)),
                0, false);

        segment.setTerminalStart(SegmentTerminal.TELOMERE);
        segment.setTerminalEnd(SegmentTerminal.TELOMERE);

        return segment;

        /*
        return ImmutableSegment.builder()
                .sampleId(sampleId)
                .clusterId(-1)
                .chainId(-1)
                .chromosome(chromosome)
                .start(1)
                .end(refGenCoords.length(chromosome))
                .track(0)
                .startTerminal(SegmentTerminal.TELOMERE)
                .endTerminal(SegmentTerminal.TELOMERE)
                .ploidy(0)
                .inDoubleMinute(false)
                .frame(0)
                .build();
        */
    }

    private static VisSegment centromere(final String sampleId, final String chromosome, final RefGenomeCoordinates refGenCoords)
    {
        int position = refGenCoords.centromere(chromosome);

        VisSegment segment = new VisSegment(
                sampleId, -1, -1, chromosome, String.valueOf(position), String.valueOf(position), 0, false);

        segment.setTerminalStart(SegmentTerminal.CENTROMERE);
        segment.setTerminalEnd(SegmentTerminal.CENTROMERE);

        return segment;

        /*
        return ImmutableSegment.builder()
                .sampleId(sampleId)
                .clusterId(-1)
                .chainId(-1)
                .chromosome(chromosome)
                .start(position)
                .end(position)
                .track(0)
                .startTerminal(SegmentTerminal.CENTROMERE)
                .endTerminal(SegmentTerminal.CENTROMERE)
                .ploidy(0)
                .inDoubleMinute(false)
                .frame(0)
                .build();
        */
    }

    public static List<VisSegment> readSegments(final String fileName) throws IOException
    {
        return VisSegment.read(fileName);
    }

    public static List<VisSegment> extendTerminals(int terminalDistance, final List<VisSegment> segments, final List<VisSvData> links,
            final List<GenomePosition> allPositions, boolean showSimpleSvSegments, final RefGenomeCoordinates refGenCoords)
    {
        final Map<String,Integer> minPositionPerChromosome = minPositionPerChromosome(allPositions);
        final Map<String,Integer> maxPositionPerChromosome = maxPositionPerChromosome(allPositions);

        final List<VisSegment> result = Lists.newArrayList();

        for(VisSegment segment : segments)
        {
            int centromere = refGenCoords.centromere(segment.chromosome());
            int length = refGenCoords.length(segment.chromosome());

            if(segment.startTerminal() != SegmentTerminal.NONE)
            {
                final int minPositionOnChromosome = minPositionPerChromosome.get(segment.chromosome());

                final int startPosition = segment.startTerminal() == SegmentTerminal.CENTROMERE && minPositionOnChromosome < centromere
                        ? centromere
                        : Math.max(1, minPositionOnChromosome - terminalDistance);

                segment = VisSegment.from(segment);
                segment.setStart(startPosition);

                // segment = ImmutableSegment.builder().from(segment).start(startPosition).build();
            }

            if(segment.endTerminal() != SegmentTerminal.NONE)
            {
                final int maxPositionOnChromosome = maxPositionPerChromosome.get(segment.chromosome());
                final int endPosition = segment.endTerminal() == SegmentTerminal.CENTROMERE && maxPositionOnChromosome > centromere
                        ? centromere
                        : Math.min(length, maxPositionOnChromosome + terminalDistance);

                segment = VisSegment.from(segment);
                segment.setEnd(endPosition);

                // segment = ImmutableSegment.builder().from(segment).end(endPosition).build();
            }

            result.add(segment);
        }

        return incrementOnChromosome(addCentromeres(result, refGenCoords), links, showSimpleSvSegments);
    }

    public static List<VisSegment> addCentromeres(final List<VisSegment> segments, final RefGenomeCoordinates refGenCoords)
    {
        if(segments.isEmpty())
            return segments;

        final List<VisSegment> result = Lists.newArrayList(segments);
        final Set<String> existingCentromeres = segments.stream()
                .filter(x -> x.startTerminal() == SegmentTerminal.CENTROMERE || x.endTerminal() == SegmentTerminal.CENTROMERE)
                .map(GenomeRegion::chromosome)
                .collect(Collectors.toSet());

        final Set<String> requiredCentromeres = Sets.newHashSet();
        final List<GenomeRegion> segmentSpan = Span.spanRegions(segments);
        for(final GenomeRegion genomeRegion : segmentSpan)
        {
            int centromere = refGenCoords.centromere(genomeRegion.chromosome());
            if(genomeRegion.start() < centromere && genomeRegion.end() > centromere)
            {
                requiredCentromeres.add(genomeRegion.chromosome());
            }
        }

        requiredCentromeres.removeAll(existingCentromeres);
        if(!requiredCentromeres.isEmpty())
        {
            final String sampleId = segments.get(0).SampleId;
            for(final String requiredCentromere : requiredCentromeres)
            {
                result.add(centromere(sampleId, requiredCentromere, refGenCoords));
            }
        }

        return result;
    }

    public static List<VisSegment> incrementOnChromosome(final List<VisSegment> segments, final List<VisSvData> links, boolean showSimpleSvSegments)
    {
        final Map<String, Integer> trackMap = Maps.newHashMap();
        final List<VisSegment> result = Lists.newArrayList();

        int frame = -1;
        int currentTrack = 1;
        for(final VisSegment segment : segments)
        {
            if(segment.ClusterId == -1)
            {
                VisSegment newSegment = VisSegment.from(segment);
                newSegment.Track = 0;
                result.add(newSegment);
            }
            else if(showSegment(showSimpleSvSegments, segment, links))
            {
                final String chromosome = segment.chromosome();
                if(!trackMap.containsKey(chromosome))
                {
                    currentTrack = 1;
                }
                else
                {
                    currentTrack = trackMap.get(chromosome) + 1;
                }
                trackMap.put(chromosome, currentTrack);

                VisSegment newSegment = VisSegment.from(segment);
                newSegment.Track = currentTrack;
                newSegment.Frame += 2;
                result.add(newSegment);
                //result.add(ImmutableSegment.builder().from(segment).track(currentTrack).frame(frame += 2).build());
            }
        }

        return result;
    }

    private static boolean showSegment(boolean showSimpleSvSegments, final VisSegment segment, final List<VisSvData> links)
    {
        if(segment.startTerminal() == SegmentTerminal.NONE && segment.endTerminal() == SegmentTerminal.NONE)
        {
            return true;
        }

        GenomePosition startPosition = GenomePositions.create(segment.chromosome(), segment.start());
        if(segment.startTerminal() == SegmentTerminal.NONE && showSegment(showSimpleSvSegments, startPosition, links))
        {
            return true;
        }

        GenomePosition endPosition = GenomePositions.create(segment.chromosome(), segment.end());
        return segment.endTerminal() == SegmentTerminal.NONE && showSegment(showSimpleSvSegments, endPosition, links);

    }

    private static boolean showSegment(boolean showSimpleSvSegments, final GenomePosition position, final List<VisSvData> links)
    {
        return VisLinks.findLink(position, links).filter(x -> !x.connectorsOnly(showSimpleSvSegments)).isPresent();
    }

}
