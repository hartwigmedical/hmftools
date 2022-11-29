package com.hartwig.hmftools.linx.visualiser.circos;

import static com.hartwig.hmftools.linx.analysis.SvUtilities.loadConfigFile;
import static com.hartwig.hmftools.linx.annotators.FragileSiteAnnotator.fragileSitesResourceFile;
import static com.hartwig.hmftools.linx.annotators.LineElementAnnotator.lineElementsResourceFile;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;

public class Highlights
{
    public static final List<GenomeRegion> FRAGILE_SITES = Lists.newArrayList();
    public static final List<GenomeRegion> LINE_ELEMENTS = Lists.newArrayList();

    public static List<GenomeRegion> limitHighlightsToRegions(final List<GenomeRegion> highlights, final List<GenomeRegion> segments)
    {
        final List<GenomeRegion> result = Lists.newArrayList();

        for (final GenomeRegion highlight : highlights)
        {
            final String contig = highlight.chromosome();
            final List<GenomeRegion> chromosomeSegments =
                    segments.stream().filter(x -> x.chromosome().equals(contig)).collect(Collectors.toList());
            if (!chromosomeSegments.isEmpty())
            {
                int minTrackPosition = chromosomeSegments.stream().mapToInt(GenomeRegion::start).min().orElse(0);
                int maxTrackPosition = chromosomeSegments.stream().mapToInt(GenomeRegion::end).max().orElse(0);
                if (highlight.end() >= minTrackPosition && highlight.start() <= maxTrackPosition)
                {

                    result.add(GenomeRegions.create(contig,
                            Math.max(minTrackPosition, highlight.start()),
                            Math.min(maxTrackPosition, highlight.end())));
                }
            }
        }
        return result;
    }

    public static void populateKnownSites(final RefGenomeVersion refGenomeVersion)
    {
        List<String> resourceLines = new BufferedReader(new InputStreamReader(
                Highlights.class.getResourceAsStream(fragileSitesResourceFile(refGenomeVersion))))
                .lines().collect(Collectors.toList());

        loadConfigFile(resourceLines, refGenomeVersion).forEach(x -> FRAGILE_SITES.add(x.genomeRegion()));

        resourceLines = new BufferedReader(new InputStreamReader(
                Highlights.class.getResourceAsStream(lineElementsResourceFile(refGenomeVersion))))
                .lines().collect(Collectors.toList());

        loadConfigFile(resourceLines, refGenomeVersion).forEach(x -> LINE_ELEMENTS.add(x.genomeRegion()));
    }
}
