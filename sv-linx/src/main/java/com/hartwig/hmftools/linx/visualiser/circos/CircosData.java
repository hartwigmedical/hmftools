package com.hartwig.hmftools.linx.visualiser.circos;

import static java.util.stream.Collectors.toList;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.linx.visualiser.data.CopyNumberAlteration;
import com.hartwig.hmftools.linx.visualiser.data.Exon;
import com.hartwig.hmftools.linx.visualiser.data.Exons;
import com.hartwig.hmftools.linx.visualiser.data.Link;
import com.hartwig.hmftools.linx.visualiser.data.Links;
import com.hartwig.hmftools.linx.visualiser.data.Segment;

import org.jetbrains.annotations.NotNull;

public class CircosData
{
    @NotNull
    private final List<Segment> unadjustedSegments;
    @NotNull
    private final List<Link> unadjustedLinks;
    @NotNull
    private final List<CopyNumberAlteration> unadjustedAlterations;
    @NotNull
    private final List<Exon> unadjustedExons;
    @NotNull
    private final List<GenomeRegion> unadjustedFragileSites;
    @NotNull
    private final List<GenomeRegion> unadjustedLineElements;

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
    private final Map<String, Integer> contigLengths;

    public CircosData(@NotNull final List<Segment> unadjustedSegments,
            @NotNull final List<Link> unadjustedLinks,
            @NotNull final List<CopyNumberAlteration> unadjustedAlterations,
            @NotNull final List<Exon> unadjustedExons)
    {
        this.unadjustedSegments = unadjustedSegments;
        this.unadjustedLinks = unadjustedLinks;
        this.unadjustedAlterations = unadjustedAlterations;
        this.unadjustedExons = unadjustedExons;

        this.unadjustedFragileSites =
                Highlights.limitHighlightsToSegments(Highlights.fragileSites(), unadjustedSegments);

        this.unadjustedLineElements =
                Highlights.limitHighlightsToSegments(Highlights.lineElements(), unadjustedSegments);

        // Note we do not add exons here because we want them interpolated.
        final List<GenomePosition> unadjustedPositions = Lists.newArrayList();
        unadjustedPositions.addAll(Links.allPositions(unadjustedLinks));
        unadjustedPositions.addAll(Span.allPositions(unadjustedSegments));
        unadjustedPositions.addAll(Span.allPositions(unadjustedAlterations));
        unadjustedPositions.addAll(Span.allPositions(unadjustedFragileSites));
        unadjustedPositions.addAll(Span.allPositions(unadjustedLineElements));
        unadjustedPositions.addAll(Span.allPositions(Exons.geneSpanPerChromosome(unadjustedExons)));
//        unadjustedPositions.addAll(Span.allPositions(unadjustedExons));

        final ScalePosition scalePosition = new ScalePosition(unadjustedPositions);
        final List<GenomePosition> scaledPositions = scalePosition.scaled();
        contigLengths = contigLengths(scaledPositions);

        segments = scalePosition.scaleSegments(unadjustedSegments);
        links = scalePosition.scaleLinks(unadjustedLinks);
        alterations = scalePosition.scaleAlterations(unadjustedAlterations);
        fragileSites = scalePosition.scaleRegions(unadjustedFragileSites);
        lineElements = scalePosition.scaleRegions(unadjustedLineElements);
        exons = scalePosition.interpolateExons(unadjustedExons);

    }

    @NotNull
    public List<Segment> unadjustedSegments()
    {
        return unadjustedSegments;
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
    public List<Exon> unadjustedExons()
    {
        return unadjustedExons;
    }

    @NotNull
    public List<GenomeRegion> unadjustedFragileSites()
    {
        return unadjustedFragileSites;
    }

    @NotNull
    public List<GenomeRegion> unadjustedLineElements()
    {
        return unadjustedLineElements;
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
