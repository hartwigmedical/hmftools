package com.hartwig.hmftools.linx.visualiser.circos;

import static java.util.stream.Collectors.toList;
import static java.util.stream.Collectors.toSet;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.linx.visualiser.data.CopyNumberAlteration;
import com.hartwig.hmftools.linx.visualiser.data.Exon;
import com.hartwig.hmftools.linx.visualiser.data.Exons;
import com.hartwig.hmftools.linx.visualiser.data.Fusion;
import com.hartwig.hmftools.linx.visualiser.data.Gene;
import com.hartwig.hmftools.linx.visualiser.data.Genes;
import com.hartwig.hmftools.linx.visualiser.data.Link;
import com.hartwig.hmftools.linx.visualiser.data.Links;
import com.hartwig.hmftools.linx.visualiser.data.ProteinDomain;
import com.hartwig.hmftools.linx.visualiser.data.ProteinDomains;
import com.hartwig.hmftools.linx.visualiser.data.Segment;

import org.jetbrains.annotations.NotNull;

public class CircosData
{
    private final List<Exon> exons;
    private final List<Link> links;
    private final List<Gene> genes;
    private final List<Segment> segments;
    private final List<GenomeRegion> lineElements;
    private final List<GenomeRegion> fragileSites;
    private final List<ProteinDomain> proteinDomains;
    private final List<CopyNumberAlteration> alterations;

    private final List<Link> unadjustedLinks;
    private final List<CopyNumberAlteration> unadjustedAlterations;

    private final Map<String, Integer> contigLengths;

    private final Set<String> upstreamGenes;
    private final Set<String> downstreamGenes;

    private final int maxTracks;
    private final double maxCopyNumber;
    private final double maxMinorAllelePloidy;

    public CircosData(boolean scaleExons,
            @NotNull final List<Segment> unadjustedSegments,
            @NotNull final List<Link> unadjustedLinks,
            @NotNull final List<CopyNumberAlteration> unadjustedAlterations,
            @NotNull final List<Exon> unadjustedExons,
            @NotNull final List<ProteinDomain> unadjustedProteinDomains,
            @NotNull final List<Fusion> fusions)
    {
        this.upstreamGenes = fusions.stream().map(Fusion::geneUp).collect(toSet());
        this.downstreamGenes = fusions.stream().map(Fusion::geneDown).collect(toSet());
        this.unadjustedLinks = unadjustedLinks;
        this.unadjustedAlterations = unadjustedAlterations;
        final List<GenomeRegion> unadjustedFragileSites =
                Highlights.limitHighlightsToSegments(Highlights.fragileSites(), unadjustedSegments);

        final List<GenomeRegion> unadjustedLineElements =
                Highlights.limitHighlightsToSegments(Highlights.lineElements(), unadjustedSegments);

        final List<Gene> unadjustedGenes = Genes.uniqueGenes(unadjustedExons);
        final List<Exon> unadjustedGeneExons = Exons.geneExons(unadjustedGenes, unadjustedExons);
        final List<ProteinDomain> unadjustedGeneProteinDomains =
                ProteinDomains.geneProteinDomains(unadjustedGenes, unadjustedProteinDomains);

        final List<GenomePosition> unadjustedPositions = Lists.newArrayList();
        unadjustedPositions.addAll(Links.allPositions(unadjustedLinks));
        unadjustedPositions.addAll(Span.allPositions(unadjustedSegments));
        unadjustedPositions.addAll(Span.allPositions(unadjustedAlterations));
        unadjustedPositions.addAll(Span.allPositions(unadjustedFragileSites));
        unadjustedPositions.addAll(Span.allPositions(unadjustedLineElements));
        unadjustedGenes.stream().map(x -> GenomePositions.create(x.chromosome(), x.namePosition())).forEach(unadjustedPositions::add);
        if (scaleExons)
        {
            unadjustedPositions.addAll(Span.allPositions(unadjustedGeneProteinDomains));
            unadjustedPositions.addAll(Span.allPositions(unadjustedGeneExons));
        }
        else
        {
            unadjustedPositions.addAll(Span.allPositions(Exons.geneSpanPerChromosome(unadjustedGeneExons)));
        }

        final ScalePosition scalePosition = new ScalePosition(unadjustedPositions);
        final List<GenomePosition> scaledPositions = scalePosition.scaled();
        contigLengths = contigLengths(scaledPositions);

        segments = scalePosition.scaleSegments(unadjustedSegments);
        links = scalePosition.scaleLinks(unadjustedLinks);
        alterations = scalePosition.scaleAlterations(unadjustedAlterations);
        fragileSites = scalePosition.scaleRegions(unadjustedFragileSites);
        lineElements = scalePosition.scaleRegions(unadjustedLineElements);
        genes = scalePosition.scaleGene(unadjustedGenes);

        // Note the following *might* be interpolated
        exons = scalePosition.interpolateExons(unadjustedGeneExons);
        proteinDomains = scalePosition.interpolateProteinDomains(unadjustedGeneProteinDomains);

        maxTracks = segments.stream().mapToInt(Segment::track).max().orElse(0) + 1;
        maxCopyNumber = alterations.stream().mapToDouble(CopyNumberAlteration::copyNumber).max().orElse(0);
        maxMinorAllelePloidy = alterations.stream().mapToDouble(CopyNumberAlteration::minorAllelePloidy).max().orElse(0);
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
        return !exons.isEmpty();
    }

    public int maxTracks()
    {
        return maxTracks;
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
    public List<Link> unadjustedLinks()
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
    public List<Link> links()
    {
        return links;
    }

    @NotNull
    public List<CopyNumberAlteration> alterations()
    {
        return alterations;
    }

    @NotNull
    public List<Exon> exons()
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
    public List<ProteinDomain> proteinDomains()
    {
        return proteinDomains;
    }

    @NotNull
    public Map<String, Integer> contigLengths()
    {
        return contigLengths;
    }

    public int totalContigLength()
    {
        return contigLengths().values().stream().mapToInt(x -> x).sum();
    }

    public long untruncatedCopyNumberAlterationsCount()
    {
        return alterations.stream().filter(x -> !x.truncated()).count();
    }

    @NotNull
    private Map<String, Integer> contigLengths(@NotNull final List<GenomePosition> positions)
    {
        final Map<String, Integer> results = new LinkedHashMap<>();
        final List<GenomePosition> sortedPositions = positions.stream().sorted().collect(toList());

        for (GenomePosition position : sortedPositions)
        {
            int end = (int) position.position();
            results.merge(position.chromosome(), end, Math::max);
        }
        return results;
    }
}
