package com.hartwig.hmftools.linx.visualiser.circos;

import static java.util.stream.Collectors.toSet;

import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.ImmutableCobaltRatio;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.purple.PurpleSegment;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.linx.visualiser.CircosConfig;
import com.hartwig.hmftools.linx.visualiser.data.Connector;
import com.hartwig.hmftools.linx.visualiser.data.Connectors;
import com.hartwig.hmftools.linx.visualiser.data.DisruptedExons;
import com.hartwig.hmftools.linx.visualiser.data.GeneUtils;
import com.hartwig.hmftools.linx.visualiser.data.VisExons;
import com.hartwig.hmftools.linx.visualiser.data.Gene;
import com.hartwig.hmftools.linx.visualiser.data.VisLinks;
import com.hartwig.hmftools.linx.visualiser.file.VisCopyNumber;
import com.hartwig.hmftools.linx.visualiser.file.VisFusion;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExon;
import com.hartwig.hmftools.linx.visualiser.file.VisSegment;
import com.hartwig.hmftools.linx.visualiser.file.VisSvData;

public class CircosData
{
    public final List<VisGeneExon> Exons;
    public final List<VisSvData> SvData;
    public final List<Gene> Genes;
    public final List<VisSegment> Segments;
    public final List<GenomeRegion> CentromereSites;
    public final List<GenomeRegion> LineElements;
    public final List<GenomeRegion> FragileSites;
    public final List<VisCopyNumber> CopyNumbers;
    public final List<GenomeRegion> DisruptedGeneRegions;
    public final List<Connector> Connectors;

    public final List<AmberBAF> AmberBAFs;
    public final List<CobaltRatio> CobaltRatios;
    public final List<PurpleSegment> PurpleSegments;

    public final List<VisSvData> UnadjustedLinks;
    public final List<VisCopyNumber> UnadjustedCopyNumbers;

    public final List<GenomeRegion> ChromosomeRanges;

    public final Set<GenomePosition> ContigLengths;

    public final Set<String> UpstreamGenes;
    public final Set<String> DownstreamGenes;

    public final CircosConfig Config;

    public static final int COPY_NUMBER_BASELINE = 0;
    public static final int COPY_NUMBER_LOSS_MIN = -2;
    public final double CopyNumberMax;
    public final int CopyNumberTracksMax;

    public static final int MINOR_ALLELE_COPY_NUMBER_LOSS_MIN = -1;
    public final double MinorAlleleCopyNumberMax;
    public final int MinorAlleleCopyNumberTracksMax;

    public final int SvTracksMax;

    public static final int AMBER_BAF_MIN = 0;
    public static final int AMBER_BAF_MAX = 1;

    public static final int COBALT_RATIO_MIN = 0;
    public static final int COBALT_RATIO_MAX = 2;

    public final double LabelSize;
    public final double GeneLabelSize;
    public final int MaxFrame;

    public CircosData(
            final CircosConfig config, final List<VisSegment> unadjustedSegments,
            final List<VisSvData> unadjustedLinks, final List<VisCopyNumber> unadjustedCopyNumbers,
            final List<VisGeneExon> unadjustedExons, final List<VisFusion> fusions,
            final List<AmberBAF> unadjustedAmberBAFs,
            final List<CobaltRatio> unadjustedCobaltRatios,
            final List<PurpleSegment> unadjustedPurpleSegments,
            boolean showSimpleSvSegments, boolean showFragileSites, boolean showLineElements
    )
    {
        UpstreamGenes = fusions.stream().map(x -> x.GeneNameUp).collect(toSet());
        DownstreamGenes = fusions.stream().map(x -> x.GeneNameDown).collect(toSet());
        UnadjustedLinks = unadjustedLinks;
        UnadjustedCopyNumbers = unadjustedCopyNumbers;
        Config = config;

        List<GenomeRegion> unadjustedDisruptedGeneRegions = DisruptedExons.disruptedGeneRegions(fusions, unadjustedExons);

        List<Gene> unadjustedGenes = GeneUtils.uniqueGenes(unadjustedExons);

        List<VisGeneExon> unadjustedGeneExons = VisExons.geneExons(unadjustedGenes, unadjustedExons);
        List<GenomeRegion> unadjustedGeneExonRegions = unadjustedGeneExons.stream().collect(Collectors.toList());

        // Determine the segment breakpoints of in the circos plot
        List<GenomePosition> positionsToScale = Lists.newArrayList();
        positionsToScale.addAll(VisLinks.allPositions(unadjustedLinks));
        positionsToScale.addAll(Span.allPositions(unadjustedSegments));
        positionsToScale.addAll(config.InterpolateCopyNumberPositions
                ? Span.minMaxPositions(unadjustedCopyNumbers)
                : Span.allPositions(unadjustedCopyNumbers));
        if(!config.InterpolateExonPositions)
        {
            positionsToScale.addAll(Span.allPositions(unadjustedGeneExonRegions));
        }

        positionsToScale = positionsToScale.stream().distinct().collect(Collectors.toList());
        Collections.sort(positionsToScale);

        // Optional annotations
        List<GenomeRegion> unadjustedCentromereSites =
                Highlights.limitHighlightsToRegions(Highlights.CENTROMERES, Span.spanPositions(positionsToScale));

        List<GenomeRegion> unadjustedFragileSites = !showFragileSites
                ? Lists.newArrayList()
                : Highlights.limitHighlightsToRegions(Highlights.FRAGILE_SITES, Span.spanPositions(positionsToScale));

        List<GenomeRegion> unadjustedLineElements = !showLineElements
                ? Lists.newArrayList()
                : Highlights.limitHighlightsToRegions(Highlights.LINE_ELEMENTS, Span.spanPositions(positionsToScale));

        ChromosomeRanges = Span.spanPositions(positionsToScale);

        // Scale positions based on the log of the size of each plot contig
        final ScalePosition positionScaler = new ScalePosition(positionsToScale);
        ContigLengths = positionScaler.contigLengths();
        Segments = positionScaler.scaleSegments(unadjustedSegments);
        SvData = positionScaler.scaleLinks(unadjustedLinks);
        CopyNumbers = positionScaler.interpolateCopyNumbers(unadjustedCopyNumbers);
        CentromereSites = positionScaler.interpolateRegions(unadjustedCentromereSites);
        LineElements = positionScaler.interpolateRegions(unadjustedLineElements);
        FragileSites = positionScaler.interpolateRegions(unadjustedFragileSites);
        Genes = positionScaler.interpolateGene(unadjustedGenes);
        DisruptedGeneRegions = positionScaler.interpolateRegions(unadjustedDisruptedGeneRegions);
        Exons = positionScaler.interpolateExons(unadjustedGeneExons);

        List<AmberBAF> amberBAFsDownsampled = Downsampler.downsampleWithMinimumPerContig(unadjustedAmberBAFs, positionsToScale);
        List<AmberBAF> amberBAFsScaled = positionScaler.interpolateAmberBAFs(amberBAFsDownsampled);
        AmberBAFs = amberBAFsScaled;

        List<CobaltRatio> cobaltRatiosDownsampled = Downsampler.downsampleWithMinimumPerContig(unadjustedCobaltRatios, positionsToScale);
        List<CobaltRatio> cobaltRatiosBreakpointAligned = alignCobaltPositionsToBreakpoints(cobaltRatiosDownsampled, positionsToScale);
        List<CobaltRatio> cobaltRatiosScaled = positionScaler.interpolateCobaltRatios(cobaltRatiosBreakpointAligned);
        CobaltRatios = cobaltRatiosScaled;

        PurpleSegments = positionScaler.interpolatePurpleSegments(unadjustedPurpleSegments);

        SvTracksMax = Segments.stream().mapToInt(x -> x.Track).max().orElse(0) + 1;
        CopyNumberMax = CopyNumbers.stream().mapToDouble(x -> x.CopyNumber).max().orElse(0);
        CopyNumberTracksMax = Math.max(2, (int) Math.round(Math.ceil(CopyNumberMax - 2)));
        MinorAlleleCopyNumberMax = CopyNumbers.stream().mapToDouble(VisCopyNumber::minorAlleleCopyNumber).max().orElse(0);
        MinorAlleleCopyNumberTracksMax = Math.max(1, (int) Math.round(Math.ceil(MinorAlleleCopyNumberMax - 1)));

        Connectors = new Connectors(showSimpleSvSegments).createConnectors(Segments, SvData);
        LabelSize = config.labelSize(untruncatedVisCopyNumberFilesCount());

        int actualMaxGeneCharacters = Genes.stream().mapToInt(x -> x.name().length()).max().orElse(0);
        GeneLabelSize = actualMaxGeneCharacters > config.MaxGeneCharacters
                ? 0.9d * config.MaxGeneCharacters / actualMaxGeneCharacters * LabelSize
                : LabelSize;

        MaxFrame = Segments.stream().mapToInt(x -> x.Frame).max().orElse(0);
    }

    private static final int COBALT_WINDOW_SIZE = 1000; // Importing CobaltConstants would lead to a dependency on Cobalt
    private static List<CobaltRatio> alignCobaltPositionsToBreakpoints(List<CobaltRatio> cobaltRatios, List<GenomePosition> positionsToScale)
    {
        /*
        Cobalt positions refer to the start position of the window. Change the Cobalt positions such that:

        If a window overlaps 2 segments, use the middle breakpoint as the Cobalt position
        ----|----|----|----
                 ^

        If a window overlaps an odd number of segments: place the Cobalt position at the midpoint of the middle segment
        ----|----|----
               ^ (midpoint is right aligned)
         */

        List<CobaltRatio> realignedCobaltRatios = Lists.newArrayList();
        for(CobaltRatio cobaltRatio : cobaltRatios)
        {
            GenomeRegion cobaltRegion = GenomeRegions.create(
                    cobaltRatio.chromosome(),
                    cobaltRatio.position(),
                    cobaltRatio.position() + COBALT_WINDOW_SIZE
            );

            List<GenomePosition> overlappingPositions = Lists.newArrayList();
            for(GenomePosition position : positionsToScale)
            {
                if(cobaltRegion.contains(position))
                    overlappingPositions.add(position);
            }

            if(overlappingPositions.isEmpty())
            {
                realignedCobaltRatios.add(cobaltRatio);
                continue;
            }

            int midIndex = overlappingPositions.size() / 2;
            int midPosition;
            if(overlappingPositions.size() % 2 != 0)
            {
                // Odd no. of breakpoints = even no. of segments --> Take the breakpoint
                midPosition = overlappingPositions.get(midIndex).position();
            }
            else
            {
                // Even no. of breakpoints = odd no. of segments --> Take midpoint of the middle segment
                int midIndexLeft = midIndex - 1;
                midPosition = (overlappingPositions.get(midIndexLeft).position() + overlappingPositions.get(midIndex).position()) / 2;
            }

            CobaltRatio realignedCobaltRatio = ImmutableCobaltRatio.builder().from(cobaltRatio)
                    .position(midPosition)
                    .build();

            realignedCobaltRatios.add(realignedCobaltRatio);
        }

        return realignedCobaltRatios;
    }

    public List<Connector> connectors()
    {
        return Connectors;
    }

    public boolean displayGenes()
    {
        return !Exons.isEmpty() && Doubles.positive(Config.GeneRelativeSize);
    }

    // public List<VisCopyNumber> unadjustedAlterations() { return unadjustedCopyNumbers; }
    public int totalContigLength()
    {
        return ContigLengths.stream().mapToInt(x -> (int) x.position()).sum();
    }

    private long untruncatedVisCopyNumberFilesCount()
    {
        return CopyNumbers.stream().filter(x -> !x.Truncated).count();
    }

    public double labelSize()
    {
        return LabelSize;
    }

    public double geneLabelSize()
    {
        return GeneLabelSize;
    }
}
