package com.hartwig.hmftools.linx.visualiser.circos;

import static java.util.stream.Collectors.toSet;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
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
import com.hartwig.hmftools.purple.region.ObservedRegion;

public class CircosData
{
    public final List<VisGeneExon> Exons;
    public final List<VisSvData> SvData;
    public final List<Gene> Genes;
    public final List<VisSegment> Segments;
    public final List<GenomeRegion> LineElements;
    public final List<GenomeRegion> FragileSites;
    public final List<VisCopyNumber> CopyNumbers;
    public final List<GenomeRegion> DisruptedGeneRegions;
    public final List<Connector> Connectors;

    public final List<AmberBAF> AmberBAFs;
    public final List<CobaltRatio> CobaltRatios;
    public final List<ObservedRegion> PurpleSegments;

    public final List<VisSvData> UnadjustedLinks;
    public final List<VisCopyNumber> UnadjustedCopyNumbers;

    public final Set<GenomePosition> ContigLengths;

    public final Set<String> UpstreamGenes;
    public final Set<String> DownstreamGenes;

    public final CircosConfig Config;

    public final int MaxTracks;
    public final double MaxPloidy;
    public final double MaxCopyNumber;
    public final double MaxMinorAlleleCopyNumber;
    public final double LabelSize;
    public final double GeneLabelSize;
    public final int MaxFrame;

    public CircosData(
            final CircosConfig config, final List<VisSegment> unadjustedSegments,
            final List<VisSvData> unadjustedLinks, final List<VisCopyNumber> unadjustedCopyNumbers,
            final List<VisGeneExon> unadjustedExons, final List<VisFusion> fusions,
            final List<AmberBAF> unadjustedAmberBAFs,
            final List<CobaltRatio> unadjustedCobaltRatios,
            final List<ObservedRegion> unadjustedPurpleSegments,
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

        List<GenomeRegion> unadjustedFragileSites = !showFragileSites
                ? Lists.newArrayList()
                : Highlights.limitHighlightsToRegions(Highlights.FRAGILE_SITES, Span.spanPositions(positionsToScale));

        List<GenomeRegion> unadjustedLineElements = !showLineElements
                ? Lists.newArrayList()
                : Highlights.limitHighlightsToRegions(Highlights.LINE_ELEMENTS, Span.spanPositions(positionsToScale));

        final ScalePosition scalePosition = new ScalePosition(positionsToScale);
        ContigLengths = scalePosition.contigLengths();
        Segments = scalePosition.scaleSegments(unadjustedSegments);
        SvData = scalePosition.scaleLinks(unadjustedLinks);
        CopyNumbers = scalePosition.interpolateCopyNumbers(unadjustedCopyNumbers);
        FragileSites = scalePosition.interpolateRegions(unadjustedFragileSites);
        LineElements = scalePosition.interpolateRegions(unadjustedLineElements);
        Genes = scalePosition.interpolateGene(unadjustedGenes);
        DisruptedGeneRegions = scalePosition.interpolateRegions(unadjustedDisruptedGeneRegions);
        Exons = scalePosition.interpolateExons(unadjustedGeneExons);

        AmberBAFs = scalePosition.interpolateAmberBAFs(unadjustedAmberBAFs);
        CobaltRatios = scalePosition.interpolateCobaltRatios(unadjustedCobaltRatios);
        PurpleSegments = scalePosition.interpolatePurpleSegments(unadjustedPurpleSegments);

        MaxTracks = Segments.stream().mapToInt(x -> x.Track).max().orElse(0) + 1;
        MaxCopyNumber = CopyNumbers.stream().mapToDouble(x -> x.CopyNumber).max().orElse(0);
        MaxMinorAlleleCopyNumber = CopyNumbers.stream().mapToDouble(VisCopyNumber::minorAlleleCopyNumber).max().orElse(0);

        double maxLinkPloidy = SvData.stream().mapToDouble(x -> x.JCN).max().orElse(0);
        double maxSegmentsPloidy = Segments.stream().mapToDouble(x -> x.LinkPloidy).max().orElse(0);

        MaxPloidy = Math.max(maxLinkPloidy, maxSegmentsPloidy);
        Connectors = new Connectors(showSimpleSvSegments).createConnectors(Segments, SvData);
        LabelSize = config.labelSize(untruncatedVisCopyNumberFilesCount());

        int actualMaxGeneCharacters = Genes.stream().mapToInt(x -> x.name().length()).max().orElse(0);
        GeneLabelSize = actualMaxGeneCharacters > config.MaxGeneCharacters
                ? 0.9d * config.MaxGeneCharacters / actualMaxGeneCharacters * LabelSize
                : LabelSize;

        MaxFrame = Segments.stream().mapToInt(x -> x.Frame).max().orElse(0);
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
