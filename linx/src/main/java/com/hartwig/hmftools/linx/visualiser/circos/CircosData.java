package com.hartwig.hmftools.linx.visualiser.circos;

import static java.util.stream.Collectors.toSet;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.linx.visualiser.CircosConfig;
import com.hartwig.hmftools.linx.visualiser.data.Connector;
import com.hartwig.hmftools.linx.visualiser.data.Connectors;
import com.hartwig.hmftools.linx.visualiser.data.DisruptedExons;
import com.hartwig.hmftools.linx.visualiser.data.VisExons;
import com.hartwig.hmftools.linx.visualiser.data.Gene;
import com.hartwig.hmftools.linx.visualiser.data.Genes;
import com.hartwig.hmftools.linx.visualiser.data.VisLinks;
import com.hartwig.hmftools.linx.visualiser.file.VisCopyNumber;
import com.hartwig.hmftools.linx.visualiser.file.VisFusion;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExon;
import com.hartwig.hmftools.linx.visualiser.file.VisSegment;
import com.hartwig.hmftools.linx.visualiser.file.VisSvData;

import org.jetbrains.annotations.NotNull;

public class CircosData
{
    private final List<VisGeneExon> exons;
    private final List<VisSvData> links;
    private final List<Gene> genes;
    private final List<VisSegment> segments;
    private final List<GenomeRegion> lineElements;
    private final List<GenomeRegion> fragileSites;
    private final List<VisCopyNumber> copyNumbers;
    private final List<GenomeRegion> disruptedGeneRegions;
    private final List<Connector> connectors;

    private final List<VisSvData> unadjustedLinks;
    private final List<VisCopyNumber> unadjustedCopyNumbers;

    private final Set<GenomePosition> contigLengths;

    private final Set<String> upstreamGenes;
    private final Set<String> downstreamGenes;

    private final CircosConfig config;

    private final int maxTracks;
    private final double maxPloidy;
    private final double maxCopyNumber;
    private final double maxMinorAllelePloidy;
    private final double labelSize;
    private final double geneLabelSize;
    private final int maxFrame;

    public CircosData(
            boolean showSimpleSvSegments, final CircosConfig config, final List<VisSegment> unadjustedSegments,
            final List<VisSvData> unadjustedLinks, final List<VisCopyNumber> unadjustedCopyNumbers,
            final List<VisGeneExon> unadjustedExons, final List<VisFusion> fusions)
    {
        this.upstreamGenes = fusions.stream().map(x -> x.GeneNameUp).collect(toSet());
        this.downstreamGenes = fusions.stream().map(x -> x.GeneNameDown).collect(toSet());
        this.unadjustedLinks = unadjustedLinks;
        this.unadjustedCopyNumbers = unadjustedCopyNumbers;
        this.config = config;

        final List<GenomeRegion> unadjustedDisruptedGeneRegions = DisruptedExons.disruptedGeneRegions(fusions, unadjustedExons);

        final List<Gene> unadjustedGenes = Genes.uniqueGenes(unadjustedExons);

        final List<VisGeneExon> unadjustedGeneExons = VisExons.geneExons(unadjustedGenes, unadjustedExons);
        final List<GenomeRegion> unadjustedGeneExonRegions = unadjustedGeneExons.stream().collect(Collectors.toList());

        final List<GenomePosition> positionsToScale = Lists.newArrayList();
        positionsToScale.addAll(VisLinks.allPositions(unadjustedLinks));
        positionsToScale.addAll(Span.allPositions(unadjustedSegments));
        positionsToScale.addAll(config.InterpolateCopyNumberPositions
                ? Span.minMaxPositions(unadjustedCopyNumbers)
                : Span.allPositions(unadjustedCopyNumbers));
        if(!config.InterpolateExonPositions)
        {
            positionsToScale.addAll(Span.allPositions(unadjustedGeneExonRegions));
        }

        final List<GenomeRegion> unadjustedFragileSites =
                Highlights.limitHighlightsToRegions(Highlights.FRAGILE_SITES, Span.spanPositions(positionsToScale));

        final List<GenomeRegion> unadjustedLineElements =
                Highlights.limitHighlightsToRegions(Highlights.LINE_ELEMENTS, Span.spanPositions(positionsToScale));

        final ScalePosition scalePosition = new ScalePosition(positionsToScale);
        contigLengths = scalePosition.contigLengths();
        segments = scalePosition.scaleSegments(unadjustedSegments);
        links = scalePosition.scaleLinks(unadjustedLinks);
        copyNumbers = scalePosition.interpolateCopyNumbers(unadjustedCopyNumbers);
        fragileSites = scalePosition.interpolateRegions(unadjustedFragileSites);
        lineElements = scalePosition.interpolateRegions(unadjustedLineElements);
        genes = scalePosition.interpolateGene(unadjustedGenes);
        disruptedGeneRegions = scalePosition.interpolateRegions(unadjustedDisruptedGeneRegions);
        exons = scalePosition.interpolateExons(unadjustedGeneExons);

        maxTracks = segments.stream().mapToInt(x -> x.Track).max().orElse(0) + 1;
        maxCopyNumber = copyNumbers.stream().mapToDouble(x -> x.CopyNumber).max().orElse(0);
        maxMinorAllelePloidy = copyNumbers.stream().mapToDouble(VisCopyNumber::minorAlleleCopyNumber).max().orElse(0);

        double maxLinkPloidy = links.stream().mapToDouble(x -> x.JCN).max().orElse(0);
        double maxSegmentsPloidy = segments.stream().mapToDouble(x -> x.LinkPloidy).max().orElse(0);

        maxPloidy = Math.max(maxLinkPloidy, maxSegmentsPloidy);
        connectors = new Connectors(showSimpleSvSegments).createConnectors(segments, links);
        labelSize = config.labelSize(untruncatedVisCopyNumberFilesCount());

        int actualMaxGeneCharacters = genes.stream().mapToInt(x -> x.name().length()).max().orElse(0);
        geneLabelSize = actualMaxGeneCharacters > config.MaxGeneCharacters
                ? 0.9d * config.MaxGeneCharacters / actualMaxGeneCharacters * labelSize
                : labelSize;

        maxFrame = segments.stream().mapToInt(x -> x.Frame).max().orElse(0);
    }

    public List<Connector> connectors()
    {
        return connectors;
    }

    @NotNull
    public List<GenomeRegion> disruptedGeneRegions()
    {
        return disruptedGeneRegions;
    }

    @NotNull
    public Set<String> upstreamGenes()
    {
        return upstreamGenes;
    }

    @NotNull
    public Set<String> downstreamGenes()
    {
        return downstreamGenes;
    }

    public boolean displayGenes()
    {
        return !exons.isEmpty() && Doubles.positive(config.GeneRelativeSize);
    }

    public int maxTracks()
    {
        return maxTracks;
    }

    public double maxPloidy()
    {
        return maxPloidy;
    }

    public int maxFrame()
    {
        return maxFrame;
    }

    public double maxCopyNumber()
    {
        return maxCopyNumber;
    }

    public double maxMinorAllelePloidy()
    {
        return maxMinorAllelePloidy;
    }

    public List<Gene> genes() { return genes; }
    public List<VisSvData> unadjustedLinks() { return unadjustedLinks; }
    public List<VisCopyNumber> unadjustedAlterations() { return unadjustedCopyNumbers; }
    public List<VisSegment> segments() { return segments; }
    public List<VisSvData> links() { return links; }
    public List<VisCopyNumber> copyNumbers() { return copyNumbers; }
    public List<VisGeneExon> exons() { return exons; }
    public List<GenomeRegion> fragileSites() { return fragileSites; }
    public List<GenomeRegion> lineElements() { return lineElements; }

    public Set<GenomePosition> contigLengths()
    {
        return contigLengths;
    }

    public int totalContigLength()
    {
        return contigLengths().stream().mapToInt(x -> (int) x.position()).sum();
    }

    private long untruncatedVisCopyNumberFilesCount()
    {
        return copyNumbers.stream().filter(x -> !x.Truncated).count();
    }

    public double labelSize()
    {
        return labelSize;
    }

    public double geneLabelSize()
    {
        return geneLabelSize;
    }
}
