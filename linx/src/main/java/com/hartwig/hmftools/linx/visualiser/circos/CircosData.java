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
import com.hartwig.hmftools.linx.visualiser.data.CopyNumberAlteration;
import com.hartwig.hmftools.linx.visualiser.data.DisruptedExons;
import com.hartwig.hmftools.linx.visualiser.data.VisExons;
import com.hartwig.hmftools.linx.visualiser.data.Fusion;
import com.hartwig.hmftools.linx.visualiser.data.Gene;
import com.hartwig.hmftools.linx.visualiser.data.Genes;
import com.hartwig.hmftools.linx.visualiser.data.VisSvData;
import com.hartwig.hmftools.linx.visualiser.data.VisLinks;
import com.hartwig.hmftools.linx.visualiser.data.Segment;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExon;

import org.jetbrains.annotations.NotNull;

public class CircosData
{
    private final List<VisGeneExon> exons;
    private final List<VisSvData> links;
    private final List<Gene> genes;
    private final List<Segment> segments;
    private final List<GenomeRegion> lineElements;
    private final List<GenomeRegion> fragileSites;
    private final List<CopyNumberAlteration> alterations;
    private final List<GenomeRegion> disruptedGeneRegions;
    private final List<Connector> connectors;

    private final List<VisSvData> unadjustedLinks;
    private final List<CopyNumberAlteration> unadjustedAlterations;

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
            boolean showSimpleSvSegments, final CircosConfig config, final List<Segment> unadjustedSegments,
            final List<VisSvData> unadjustedLinks, final List<CopyNumberAlteration> unadjustedAlterations,
            final List<VisGeneExon> unadjustedExons, final List<Fusion> fusions)
    {
        this.upstreamGenes = fusions.stream().map(Fusion::geneUp).collect(toSet());
        this.downstreamGenes = fusions.stream().map(Fusion::geneDown).collect(toSet());
        this.unadjustedLinks = unadjustedLinks;
        this.unadjustedAlterations = unadjustedAlterations;
        this.config = config;

        final List<GenomeRegion> unadjustedDisruptedGeneRegions = DisruptedExons.disruptedGeneRegions(fusions, unadjustedExons);

        final List<Gene> unadjustedGenes = Genes.uniqueGenes(unadjustedExons);

        final List<VisGeneExon> unadjustedGeneExons = VisExons.geneExons(unadjustedGenes, unadjustedExons);
        final List<GenomeRegion> unadjustedGeneExonRegions = unadjustedGeneExons.stream().collect(Collectors.toList());

        final List<GenomePosition> positionsToScale = Lists.newArrayList();
        positionsToScale.addAll(VisLinks.allPositions(unadjustedLinks));
        positionsToScale.addAll(Span.allPositions(unadjustedSegments));
        positionsToScale.addAll(config.InterpolateCopyNumberPositions
                ? Span.minMaxPositions(unadjustedAlterations)
                : Span.allPositions(unadjustedAlterations));
        if (!config.InterpolateExonPositions)
        {
            positionsToScale.addAll(Span.allPositions(unadjustedGeneExonRegions));
        }

        final List<GenomeRegion> unadjustedFragileSites =
                Highlights.limitHighlightsToRegions(Highlights.fragileSites(), Span.spanPositions(positionsToScale));

        final List<GenomeRegion> unadjustedLineElements =
                Highlights.limitHighlightsToRegions(Highlights.lineElements(), Span.spanPositions(positionsToScale));

        final ScalePosition scalePosition = new ScalePosition(positionsToScale);
        contigLengths = scalePosition.contigLengths();
        segments = scalePosition.scaleSegments(unadjustedSegments);
        links = scalePosition.scaleLinks(unadjustedLinks);
        alterations = scalePosition.interpolateAlterations(unadjustedAlterations);
        fragileSites = scalePosition.interpolateRegions(unadjustedFragileSites);
        lineElements = scalePosition.interpolateRegions(unadjustedLineElements);
        genes = scalePosition.interpolateGene(unadjustedGenes);
        disruptedGeneRegions = scalePosition.interpolateRegions(unadjustedDisruptedGeneRegions);
        exons = scalePosition.interpolateExons(unadjustedGeneExons);

        maxTracks = segments.stream().mapToInt(Segment::track).max().orElse(0) + 1;
        maxCopyNumber = alterations.stream().mapToDouble(CopyNumberAlteration::copyNumber).max().orElse(0);
        maxMinorAllelePloidy = alterations.stream().mapToDouble(CopyNumberAlteration::minorAlleleCopyNumber).max().orElse(0);

        double maxLinkPloidy = links.stream().mapToDouble(VisSvData::jcn).max().orElse(0);
        double maxSegmentsPloidy = segments.stream().mapToDouble(Segment::ploidy).max().orElse(0);

        maxPloidy = Math.max(maxLinkPloidy, maxSegmentsPloidy);
        connectors = new Connectors(showSimpleSvSegments).createConnectors(segments, links);
        labelSize = config.labelSize(untruncatedCopyNumberAlterationsCount());

        int actualMaxGeneCharacters = genes.stream().mapToInt(x -> x.name().length()).max().orElse(0);
        geneLabelSize = actualMaxGeneCharacters > config.MaxGeneCharacters
                ? 0.9d * config.MaxGeneCharacters / actualMaxGeneCharacters * labelSize
                : labelSize;

        maxFrame = segments.stream().mapToInt(Segment::frame).max().orElse(0);
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

    @NotNull
    public List<Gene> genes()
    {
        return genes;
    }

    @NotNull
    public List<VisSvData> unadjustedLinks()
    {
        return unadjustedLinks;
    }

    @NotNull
    public List<CopyNumberAlteration> unadjustedAlterations()
    {
        return unadjustedAlterations;
    }

    @NotNull
    public List<Segment> segments()
    {
        return segments;
    }

    @NotNull
    public List<VisSvData> links()
    {
        return links;
    }

    @NotNull
    public List<CopyNumberAlteration> alterations()
    {
        return alterations;
    }

    @NotNull
    public List<VisGeneExon> exons()
    {
        return exons;
    }

    @NotNull
    public List<GenomeRegion> fragileSites()
    {
        return fragileSites;
    }

    @NotNull
    public List<GenomeRegion> lineElements()
    {
        return lineElements;
    }

    @NotNull
    public Set<GenomePosition> contigLengths()
    {
        return contigLengths;
    }

    public int totalContigLength()
    {
        return contigLengths().stream().mapToInt(x -> (int) x.position()).sum();
    }

    private long untruncatedCopyNumberAlterationsCount()
    {
        return alterations.stream().filter(x -> !x.truncated()).count();
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
