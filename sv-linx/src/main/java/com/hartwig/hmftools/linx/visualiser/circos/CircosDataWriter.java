package com.hartwig.hmftools.linx.visualiser.circos;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.linx.visualiser.SvCircosConfig;
import com.hartwig.hmftools.linx.visualiser.data.AdjustedPosition;
import com.hartwig.hmftools.linx.visualiser.data.AdjustedPositions;
import com.hartwig.hmftools.linx.visualiser.data.CopyNumberAlteration;
import com.hartwig.hmftools.linx.visualiser.data.Exon;
import com.hartwig.hmftools.linx.visualiser.data.Gene;
import com.hartwig.hmftools.linx.visualiser.data.Link;
import com.hartwig.hmftools.linx.visualiser.data.Links;
import com.hartwig.hmftools.linx.visualiser.data.ProteinDomain;
import com.hartwig.hmftools.linx.visualiser.data.Segment;

import org.jetbrains.annotations.NotNull;

public class CircosDataWriter
{
    private static final String SINGLE_BLUE = "(107,174,214)";
    private static final String SINGLE_RED = "(214,144,107)";
    private static final String SINGLE_GREEN = "(107,214,148)";

    private static final int MAX_CONTIG_LENGTH_TO_DISPLAY_EXON_RANK = 100000;

    private static final int MIN_KAROTYPE_LENGTH = 10;
    private static final String DELIMITER = "\t";

    private final String filePrefix;
    private final ColorPicker colorPicker;
    private final SvCircosConfig circosConfig;
    private final CircosConfigWriter configWriter;
    private final ProteinDomainColors proteinDomainColors;

    public CircosDataWriter(final ColorPicker colorPicker, @NotNull final String sample,
            @NotNull final String outputDir,
            @NotNull final SvCircosConfig circosConfig,
            @NotNull final CircosConfigWriter configWriter,
            @NotNull final ProteinDomainColors proteinDomainColors)
    {
        this.colorPicker = colorPicker;
        this.configWriter = configWriter;
        this.circosConfig = circosConfig;
        this.proteinDomainColors = proteinDomainColors;
        this.filePrefix = outputDir + File.separator + sample;
    }

    public void write(@NotNull final CircosData data) throws IOException
    {
        final Map<String, String> geneColorMap = Maps.newHashMap();
        data.genes().forEach(x -> geneColorMap.put(x.name(), SINGLE_GREEN));
        data.upstreamGenes().forEach(x -> geneColorMap.put(x, SINGLE_BLUE));
        data.downstreamGenes().forEach(x -> geneColorMap.put(x, SINGLE_RED));

        final Map<String, Integer> contigLengths = data.contigLengths();
        int totalContigLength = data.totalContigLength();

        final List<Segment> segments = data.segments();
        final List<Link> links = data.links();
        final List<CopyNumberAlteration> alterations = data.alterations();
        final List<GenomeRegion> fragileSites = data.fragileSites();
        final List<GenomeRegion> lineElements = data.lineElements();
        final List<Exon> exons = data.exons();

        final String proteinDomainPath = filePrefix + ".protein_domain.circos";
        Files.write(new File(proteinDomainPath).toPath(), proteinDomain(data.proteinDomains()));

        final String exonPath = filePrefix + ".exon.circos";
        Files.write(new File(exonPath).toPath(), exons(geneColorMap, exons));

        final String exonRankPath = filePrefix + ".exon.rank.circos";
        Files.write(new File(exonRankPath).toPath(), exonRank(totalContigLength, exons));

        final String genePath = filePrefix + ".gene.circos";
        Files.write(new File(genePath).toPath(), genes(geneColorMap, data.genes()));

        final String geneNamePath = filePrefix + ".gene.name.circos";
        Files.write(new File(geneNamePath).toPath(), geneName(data.genes()));

        final String textPath = filePrefix + ".position.circos";
        Files.write(new File(textPath).toPath(), createPositionText(data.unadjustedLinks(), links));

        final String histogramPath = filePrefix + ".segment.circos";
        Files.write(new File(histogramPath).toPath(), createHistogramTrack(segments));

        final String karyotypePath = filePrefix + ".karyotype.circos";
        Files.write(new File(karyotypePath).toPath(), createKaryotypes(contigLengths));

        final String connectorPath = filePrefix + ".connector.circos";
        Files.write(new File(connectorPath).toPath(), createConnectors(segments, links));

        final String linkPath = filePrefix + ".link.circos";
        Files.write(new File(linkPath).toPath(), createLinks(links));

        final String scatterPath = filePrefix + ".scatter.circos";
        Files.write(new File(scatterPath).toPath(), createScatter(segments, links));

        final String cnaPath = filePrefix + ".cna.circos";
        Files.write(new File(cnaPath).toPath(), createCNA(alterations));

        final String mapPath = filePrefix + ".map.circos";
        Files.write(new File(mapPath).toPath(), createMinorAllelePloidy(alterations));

        final String fragile = filePrefix + ".fragile.circos";
        Files.write(new File(fragile).toPath(), highlights(fragileSites));

        final String line = filePrefix + ".line_element.circos";
        Files.write(new File(line).toPath(), highlights(lineElements));

        final String distances = filePrefix + ".distance.circos";
        Files.write(new File(distances).toPath(), createDistances(data.unadjustedAlterations(), alterations));

        //        final String chromosomeBandPath = filePrefix + ".chromosome.circos";
        //        Files.write(new File(chromosomeBandPath).toPath(), chromosomeLocations(data.unadjustedAlterations()));

    }

    //    @NotNull
    //    private List<String> chromosomeLocations(List<CopyNumberAlteration> unadjustedAlterations) {
    //        List<String> result = Lists.newArrayList();
    //
    //        List<GenomeRegion> regions = Span.spanRegions(unadjustedAlterations);
    //        for (GenomeRegion region : regions)
    //        {
    //            final String bandString = new StringJoiner(DELIMITER).add(region.chromosome())
    //                    .add(String.valueOf(region.start()))
    //                    .add(String.valueOf(region.end()))
    //                    .toString();
    //            result.add(bandString);
    //        }
    //
    //        return result;
    //    }

    @NotNull
    private List<String> genes(@NotNull final Map<String, String> geneColours, @NotNull final List<Gene> genes)
    {
        final List<String> result = Lists.newArrayList();
        for (final Gene gene : genes)
        {
            final String exonString = new StringJoiner(DELIMITER).add(circosContig(gene.chromosome()))
                    .add(String.valueOf(gene.start()))
                    .add(String.valueOf(gene.end()))
                    .add(String.valueOf(1))
                    .add("fill_color=" + geneColours.get(gene.name()))
                    .toString();
            result.add(exonString);

        }

        return result;
    }

    @NotNull
    private List<String> proteinDomain(@NotNull final List<ProteinDomain> proteinDomains)
    {
        final List<String> result = Lists.newArrayList();
        for (final ProteinDomain proteinDomain : proteinDomains)
        {
            final String color = proteinDomainColors.rgb(proteinDomain.name());

            final String exonString = new StringJoiner(DELIMITER).add(circosContig(proteinDomain.chromosome()))
                    .add(String.valueOf(proteinDomain.start()))
                    .add(String.valueOf(proteinDomain.end()))
                    .add(String.valueOf(1))
                    .add("fill_color=" + color + ",name=" + proteinDomain.name().replace(' ', '.'))
                    .toString();
            result.add(exonString);

        }

        return result;
    }

    @NotNull
    private List<String> exonRank(int totalContigLength, @NotNull final List<Exon> exons)
    {
        final List<String> result = Lists.newArrayList();
        if (totalContigLength <= MAX_CONTIG_LENGTH_TO_DISPLAY_EXON_RANK)
        {
            for (final Exon exon : exons)
            {
                long position = exon.start() + (exon.end() - exon.start()) / 2;

                final String exonString = new StringJoiner(DELIMITER).add(circosContig(exon.chromosome()))
                        .add(String.valueOf(position))
                        .add(String.valueOf(position))
                        .add(String.valueOf(exon.rank()))
                        .toString();
                result.add(exonString);

            }
        }

        return result;
    }

    @NotNull
    private List<String> geneName(@NotNull final List<Gene> genes)
    {
        final List<String> result = Lists.newArrayList();
        for (final Gene gene : genes)
        {
            final double labelSize = geneNameLabelSize(gene.name());

            final String exonString = new StringJoiner(DELIMITER).add(circosContig(gene.chromosome()))
                    .add(String.valueOf(gene.namePosition()))
                    .add(String.valueOf(gene.namePosition()))
                    .add(gene.name())
                    .toString();
            result.add(exonString);
        }

        return result;
    }

    private static double geneNameLabelSize(@NotNull final String gene)
    {
        double availablePixels = CircosConfigWriter.PIXELS * (1 - CircosConfigWriter.EXON_INNER_RADIUS);
        return Math.min(26, 4 + Math.floor(availablePixels / gene.length()));
    }

    @NotNull
    private List<String> exons(@NotNull final Map<String, String> geneColours, @NotNull final List<Exon> exons)
    {
        final List<String> result = Lists.newArrayList();
        for (final Exon exon : exons)
        {
            final String exonString = new StringJoiner(DELIMITER).add(circosContig(exon.chromosome()))
                    .add(String.valueOf(exon.start()))
                    .add(String.valueOf(exon.end()))
                    .add(String.valueOf(1))
                    .add("fill_color=" + geneColours.get(exon.gene()))
                    .toString();
            result.add(exonString);
        }

        return result;
    }

    @NotNull
    private List<String> highlights(@NotNull final List<GenomeRegion> regions)
    {
        return regions.stream()
                .map(x -> new StringJoiner(DELIMITER).add(circosContig(x.chromosome()))
                        .add(String.valueOf(x.start()))
                        .add(String.valueOf(x.end()))
                        .toString())
                .collect(toList());
    }

    @NotNull
    private List<String> createDistances(@NotNull final List<CopyNumberAlteration> unadjustedSegment,
            @NotNull final List<CopyNumberAlteration> segments)
    {
        final List<String> result = Lists.newArrayList();
        long unadjustedSegments = segments.stream().filter(x -> !x.truncated()).count();
        if (unadjustedSegments <= circosConfig.maxDistanceLabels())
        {
            for (int i = 0; i < unadjustedSegment.size(); i++)
            {
                final CopyNumberAlteration adjusted = segments.get(i);
                final CopyNumberAlteration unadjusted = unadjustedSegment.get(i);
                if (!adjusted.truncated())
                {
                    final String distance = new StringJoiner(DELIMITER).add(circosContig(adjusted.chromosome()))
                            .add(String.valueOf(adjusted.start()))
                            .add(String.valueOf(adjusted.end()))
                            .add(shorthand(unadjusted.end() - unadjusted.start()))
                            .toString();
                    result.add(distance);
                }
            }

        }

        return result;
    }

    @NotNull
    private List<String> createCNA(@NotNull final List<CopyNumberAlteration> alterations)
    {
        final List<String> result = Lists.newArrayList();
        for (CopyNumberAlteration alteration : alterations)
        {
            final String cna = new StringJoiner(DELIMITER).add(circosContig(alteration.chromosome()))
                    .add(String.valueOf(alteration.start()))
                    .add(String.valueOf(alteration.end()))
                    .add(String.valueOf(alteration.copyNumber() - 2))
                    .toString();
            result.add(cna);
        }

        return result;
    }

    @NotNull
    private List<String> createMinorAllelePloidy(@NotNull final List<CopyNumberAlteration> alterations)
    {
        final List<String> result = Lists.newArrayList();
        for (CopyNumberAlteration alteration : alterations)
        {
            final String cna = new StringJoiner(DELIMITER).add(circosContig(alteration.chromosome()))
                    .add(String.valueOf(alteration.start()))
                    .add(String.valueOf(alteration.end()))
                    .add(String.valueOf(alteration.minorAllelePloidy() - 1))
                    .toString();
            result.add(cna);
        }

        return result;
    }

    @NotNull
    private List<String> createScatter(@NotNull final List<Segment> segments, @NotNull final List<Link> links)
    {
        final List<String> result = Lists.newArrayList();
        for (Segment segment : segments)
        {

            if (segment.track() == 0)
            {
                continue;
            }

            final String colorOption = colorPicker.color(segment.clusterId(), segment.chainId());
            final String startGlyph = scatterGlyph(true, segment, links);
            result.add(scatterEntry(true, segment, colorOption, startGlyph, 20));
            if (segment.startTerminal() == SegmentTerminal.CENTROMERE)
            {
                result.add(scatterEntry(true, segment, "color=white", startGlyph, 14));
            }

            final String endGlyph = scatterGlyph(false, segment, links);
            result.add(scatterEntry(false, segment, colorOption, endGlyph, 20));
            if (segment.endTerminal() == SegmentTerminal.CENTROMERE)
            {
                result.add(scatterEntry(false, segment, "color=white", endGlyph, 14));
            }
        }

        return result;
    }

    @NotNull
    private String scatterGlyph(boolean isStart, @NotNull final Segment segment, @NotNull final List<Link> links)
    {
        long location = isStart ? segment.start() : segment.end();
        final SegmentTerminal terminal = isStart ? segment.startTerminal() : segment.endTerminal();
        if (terminal != SegmentTerminal.NONE)
        {
            return "square";
        }

        final GenomePosition startPosition = GenomePositions.create(segment.chromosome(), location);
        final boolean isFoldback =
                Links.findStartLink(startPosition, links).filter(x -> x.startInfo().equals("FOLDBACK")).isPresent() || Links.findEndLink(
                        startPosition,
                        links).filter(x -> x.endInfo().equals("FOLDBACK")).isPresent();

        return isFoldback ? "triangle" : "circle";
    }

    @NotNull
    private String scatterEntry(boolean isStart, @NotNull final Segment segment, @NotNull final String color, @NotNull final String glyph,
            int glyph_size)
    {

        long location = isStart ? segment.start() : segment.end();

        return new StringJoiner(DELIMITER).add(circosContig(segment.chromosome()))
                .add(String.valueOf(location))
                .add(String.valueOf(location))
                .add(String.valueOf(segment.track()))
                .add(color + "," + "glyph=" + glyph + ",glyph_size=" + glyph_size)
                .toString();
    }

    @NotNull
    private List<String> createLinks(@NotNull final List<Link> links)
    {
        final List<String> result = Lists.newArrayList();
        for (final Link link : Links.clean(links))
        {

            final String linkString = new StringJoiner(DELIMITER).add(circosContig(link.startChromosome()))
                    .add(String.valueOf(link.startPosition()))
                    .add(String.valueOf(link.startPosition()))
                    .add(circosContig(link.endChromosome()))
                    .add(String.valueOf(link.endPosition()))
                    .add(String.valueOf(link.endPosition()))
                    .add(colorPicker.transparentColor(link.clusterId(), link.chainId()) + "," + thicknessString(link.ploidy()))
                    .toString();
            result.add(linkString);
        }

        return result;
    }

    @NotNull
    private List<String> createConnectors(@NotNull final List<Segment> segments, @NotNull final List<Link> links)
    {
        final List<String> result = Lists.newArrayList();

        for (Segment segment : segments)
        {

            final GenomePosition startPosition = GenomePositions.create(segment.chromosome(), segment.start());

            final double r1 = configWriter.svTrackRelative(segment.track());
            double startLinkUsage = Links.linkTraverseCount(startPosition, links);

            if (startLinkUsage > 0)
            {
                long segmentsBelow = segments.stream()
                        .filter(x -> x.chromosome().equals(segment.chromosome()) && x.start() == segment.start()
                                && x.track() < segment.track())
                        .count();

                final String start = new StringJoiner(DELIMITER).add(circosContig(segment.chromosome()))
                        .add(String.valueOf(segment.start()))
                        .add(String.valueOf(segment.start()))
                        .add("r1=" + r1 + "r," + colorPicker.transparentColor(segment.clusterId(), segment.chainId()) + ","
                                + thicknessString(startLinkUsage - segmentsBelow))
                        .toString();
                result.add(start);
            }

            final GenomePosition endPosition = GenomePositions.create(segment.chromosome(), segment.end());

            double endLinkUsage = Links.linkTraverseCount(endPosition, links);

            if (endLinkUsage > 0)
            {
                long segmentsBelow = segments.stream()
                        .filter(x -> x.chromosome().equals(segment.chromosome()) && x.end() == segment.end() && x.track() < segment.track())
                        .count();
                final String end = new StringJoiner(DELIMITER).add(circosContig(segment.chromosome()))
                        .add(String.valueOf(segment.end()))
                        .add(String.valueOf(segment.end()))
                        .add("r1=" + r1 + "r," + colorPicker.transparentColor(segment.clusterId(), segment.chainId()) + ","
                                + thicknessString(endLinkUsage - segmentsBelow))
                        .toString();
                result.add(end);
            }

        }

        double rTrack1 = configWriter.svTrackRelative(0);
        for (Link link : links)
        {
            if (link.connectorsOnly())
            {
                if (link.isValidStart())
                {
                    final String start = new StringJoiner(DELIMITER).add(circosContig(link.startChromosome()))
                            .add(String.valueOf(link.startPosition()))
                            .add(String.valueOf(link.startPosition()))
                            .add("r1=" + rTrack1 + "r," + colorPicker.transparentColor(link.clusterId(), link.chainId()) + ","
                                    + thicknessString(link.ploidy()))
                            .toString();
                    result.add(start);
                }

                if (link.isValidEnd())
                {
                    final String end = new StringJoiner(DELIMITER).add(circosContig(link.endChromosome()))
                            .add(String.valueOf(link.endPosition()))
                            .add(String.valueOf(link.endPosition()))
                            .add("r1=" + rTrack1 + "r," + colorPicker.transparentColor(link.clusterId(), link.chainId()) + ","
                                    + thicknessString(link.ploidy()))
                            .toString();
                    result.add(end);
                }

            }

        }

        return result;
    }

    @NotNull
    private List<String> createKaryotypes(@NotNull final Map<String, Integer> contigLengths)
    {
        final List<String> result = Lists.newArrayList();
        for (String contig : contigLengths.keySet())
        {

            final String start = new StringJoiner(" ").add("chr -")
                    .add(circosContig(contig))
                    .add(HumanChromosome.fromString(contig).toString())
                    .add(String.valueOf(1))
                    .add(String.valueOf(Math.max(MIN_KAROTYPE_LENGTH, contigLengths.get(contig))))
                    .add("chr" + HumanChromosome.fromString(contig).toString())
                    .toString();
            result.add(start);
        }

        return result;
    }

    @NotNull
    private List<String> createHistogramTrack(@NotNull final List<Segment> segments)
    {
        final List<String> result = Lists.newArrayList();
        for (final Segment segment : segments)
        {

            if (segment.track() > 0)
            {

                double r0 = configWriter.svTrackRelative(segment.track());
                String r0String = "r0=" + r0 + "r";
                double thickness = thicknessPixels(segment.ploidy());
                String r1String = "r1=" + r0 + "r+" + thickness + "p";

                final String entry = new StringJoiner(DELIMITER).add(circosContig(segment.chromosome()))
                        .add(String.valueOf(segment.start()))
                        .add(String.valueOf(segment.end()))
                        .add(String.valueOf(segment.track()))
                        .add("fill_" + colorPicker.color(segment.clusterId(), segment.chainId()) + "," + r0String + "," + r1String)
                        .toString();
                result.add(entry);

            }
        }

        return result;
    }

    @NotNull
    private List<String> createPositionText(@NotNull final List<Link> originalLinks, @NotNull final List<Link> scaledLinks)
    {

        if (!circosConfig.displayPosition())
        {
            return Collections.emptyList();
        }

        final Set<String> result = Sets.newHashSet();
        final List<AdjustedPosition> positions = AdjustedPositions.create(originalLinks, scaledLinks);
        final Set<String> contigs = positions.stream().map(GenomePosition::chromosome).collect(Collectors.toSet());

        for (final String contig : contigs)
        {
            long currentPosition = 0;
            for (final AdjustedPosition position : positions)
            {
                if (position.chromosome().equals(contig))
                {
                    long roundedPosition = position.unadjustedPosition() / 100_000;
                    if (roundedPosition > currentPosition)
                    {
                        final String start = new StringJoiner(DELIMITER).add(circosContig(contig))
                                .add(String.valueOf(position.position()))
                                .add(String.valueOf(position.position()))
                                .add(String.valueOf(roundedPosition / 10d + "m"))
                                .toString();

                        result.add(start);
                        currentPosition = roundedPosition;
                    }
                }
            }
        }

        return result.stream().sorted().distinct().collect(toList());
    }

    @NotNull
    private static String circosContig(@NotNull final String chromosome)
    {
        return "hs" + HumanChromosome.fromString(chromosome);
    }

    @NotNull
    private static String thicknessString(double usage)
    {
        return "thickness=" + thicknessPixels(usage);
    }

    static double thicknessPixels(double usage)
    {
        return Math.min(12, Math.max(1, Math.pow(2 * usage, 1)));
    }

    @NotNull
    static String shorthand(long value)
    {
        if (value < 100)
        {
            return String.valueOf(value);
        }

        if (value < 99_950)
        {
            return String.format("%.1fk", value / 1_000d);
        }

        return String.format("%.1fm", value / 1_000_000d);
    }

}
