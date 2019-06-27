package com.hartwig.hmftools.linx.visualiser.circos;

import static java.util.stream.Collectors.toList;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.linx.visualiser.data.CopyNumberAlteration;
import com.hartwig.hmftools.linx.visualiser.data.Exon;
import com.hartwig.hmftools.linx.visualiser.data.Exons;
import com.hartwig.hmftools.linx.visualiser.data.Gene;
import com.hartwig.hmftools.linx.visualiser.data.Genes;
import com.hartwig.hmftools.linx.visualiser.data.Link;
import com.hartwig.hmftools.linx.visualiser.data.Links;
import com.hartwig.hmftools.linx.visualiser.data.ProteinDomain;
import com.hartwig.hmftools.linx.visualiser.data.Segment;

import org.jetbrains.annotations.NotNull;

public class CircosData
{
    @NotNull
    private final List<Link> unadjustedLinks;
    @NotNull
    private final List<CopyNumberAlteration> unadjustedAlterations;

    @NotNull
    private final List<Segment> segments;
    @NotNull
    private final List<Link> links;
    @NotNull
    private final List<CopyNumberAlteration> alterations;
    @NotNull
    private final List<Exon> exons;
    @NotNull
    private final List<GenomeRegion> fragileSites;
    @NotNull
    private final List<GenomeRegion> lineElements;
    @NotNull
    private final List<Gene> genes;
    @NotNull
    private  final List<ProteinDomain> proteinDomains;

    @NotNull
    private final Map<String, Integer> contigLengths;

    public CircosData(boolean scaleExons,
            @NotNull final List<Segment> unadjustedSegments,
            @NotNull final List<Link> unadjustedLinks,
            @NotNull final List<CopyNumberAlteration> unadjustedAlterations,
            @NotNull final List<Exon> unadjustedExons,
            @NotNull final List<ProteinDomain> unadjustedProteinDomains)
    {
        this.unadjustedLinks = unadjustedLinks;
        this.unadjustedAlterations = unadjustedAlterations;
        final List<GenomeRegion> unadjustedFragileSites =
                Highlights.limitHighlightsToSegments(Highlights.fragileSites(), unadjustedSegments);

        final List<GenomeRegion> unadjustedLineElements =
                Highlights.limitHighlightsToSegments(Highlights.lineElements(), unadjustedSegments);

        final List<Gene> unadjustedGenes = Genes.genes(unadjustedExons);

        final List<GenomePosition> unadjustedPositions = Lists.newArrayList();
        unadjustedPositions.addAll(Links.allPositions(unadjustedLinks));
        unadjustedPositions.addAll(Span.allPositions(unadjustedSegments));
        unadjustedPositions.addAll(Span.allPositions(unadjustedAlterations));
        unadjustedPositions.addAll(Span.allPositions(unadjustedFragileSites));
        unadjustedPositions.addAll(Span.allPositions(unadjustedLineElements));
        if (scaleExons)
        {
            unadjustedPositions.addAll(Span.allPositions(unadjustedProteinDomains));
            unadjustedPositions.addAll(Span.allPositions(unadjustedExons));
            unadjustedGenes.stream().map(x -> GenomePositions.create(x.chromosome(), x.namePosition())).forEach(unadjustedPositions::add);
        }
        else
        {
            unadjustedPositions.addAll(Span.allPositions(Exons.geneSpanPerChromosome(unadjustedExons)));
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
        exons = scalePosition.interpolateExons(unadjustedExons);
        proteinDomains = scalePosition.interpolateProteinDomains(unadjustedProteinDomains);

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
