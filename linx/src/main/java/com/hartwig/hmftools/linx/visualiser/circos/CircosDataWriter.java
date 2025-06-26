package com.hartwig.hmftools.linx.visualiser.circos;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Collections;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.purple.PurpleSegment;
import com.hartwig.hmftools.linx.visualiser.CircosConfig;
import com.hartwig.hmftools.linx.visualiser.data.AdjustedPosition;
import com.hartwig.hmftools.linx.visualiser.data.AdjustedPositions;
import com.hartwig.hmftools.linx.visualiser.data.Connector;
import com.hartwig.hmftools.linx.visualiser.data.Gene;
import com.hartwig.hmftools.linx.visualiser.file.VisCopyNumber;
import com.hartwig.hmftools.linx.visualiser.file.VisSegment;
import com.hartwig.hmftools.linx.visualiser.file.VisSvData;
import com.hartwig.hmftools.linx.visualiser.data.VisLinks;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExon;

import org.jetbrains.annotations.NotNull;

public class CircosDataWriter
{
    private static final DecimalFormat RATIO_FORMAT = new DecimalFormat("#.###", DecimalFormatSymbols.getInstance(Locale.ENGLISH));
    private static final DecimalFormat POSITION_FORMAT = new DecimalFormat("#,###");
    private static final String SINGLE_BLUE = "(107,174,214)";
    private static final String SINGLE_RED = "(214,144,107)";
    private static final String SINGLE_GREEN = "(107,214,148)";
    private static final String SINGLE_YELLOW = "(214,210,107)";

    private static final int MAX_CONTIG_LENGTH_TO_DISPLAY_EXON_RANK = 100000;

    private static final int MIN_KARYOTYPE_LENGTH = 10;
    private static final String DELIMITER = "\t";

    private final String filePrefix;
    private final ColorPicker colorPicker;
    private final CircosConfig circosConfig;
    private final CircosConfigWriter configWriter;
    private final CircosData data;
    private final Thickness thickness;

    public CircosDataWriter(final ColorPicker colorPicker, final String sample, final String outputDir,
            final CircosConfig circosConfig, final CircosConfigWriter configWriter, final CircosData data)
    {
        this.data = data;
        this.colorPicker = colorPicker;
        this.configWriter = configWriter;
        this.circosConfig = circosConfig;
        this.filePrefix = outputDir + File.separator + sample;
        this.thickness = new Thickness(circosConfig.MinLineSize, circosConfig.MaxLineSize, data.connectors());
    }

    public Object write() throws IOException
    {
        final Map<String, String> geneColorMap = Maps.newHashMap();
        data.Genes.forEach(x ->
        {
            switch (x.type())
            {
                case PSEUDOGENE:
                    geneColorMap.put(x.name(), SINGLE_YELLOW);
                    break;
                case DRIVER:
                case DISRUPTION:
                    geneColorMap.put(x.name(), SINGLE_GREEN);
                    break;
                case FUSION:
                    geneColorMap.put(x.name(), SINGLE_BLUE);
                    break;
            }
        });

        data.UpstreamGenes.forEach(x -> geneColorMap.put(x, SINGLE_BLUE));
        data.DownstreamGenes.forEach(x -> geneColorMap.put(x, SINGLE_RED));

        int totalContigLength = data.totalContigLength();

        List<VisSegment> segments = data.Segments;
        List<VisSvData> links = data.SvData;
        List<VisCopyNumber> alterations = data.CopyNumbers;
        List<GenomeRegion> centromereSites = data.CentromereSites;
        List<GenomeRegion> fragileSites = data.FragileSites;
        List<GenomeRegion> lineElements = data.LineElements;
        List<VisGeneExon> exons = data.Exons;

        String exonPath = filePrefix + ".exon.circos";
        Files.write(new File(exonPath).toPath(), exons(geneColorMap, data.DisruptedGeneRegions, exons));

        String exonRankPath = filePrefix + ".exon.rank.circos";
        Files.write(new File(exonRankPath).toPath(), exonRank(totalContigLength, exons));

        String genePath = filePrefix + ".gene.circos";
        Files.write(new File(genePath).toPath(), genes(geneColorMap, data.Genes));

        String geneNamePath = filePrefix + ".gene.name.circos";
        Files.write(new File(geneNamePath).toPath(), geneName(data.Genes));

        String textPath = filePrefix + ".position.circos";
        Files.write(new File(textPath).toPath(), createPositionText(data.UnadjustedLinks, links));

        String histogramPath = filePrefix + ".segment.circos";
        Files.write(new File(histogramPath).toPath(), createHistogramTrack(segments));

        String karyotypePath = filePrefix + ".karyotype.circos";
        Files.write(new File(karyotypePath).toPath(), createKaryotypes(data.ContigLengths));

        String connectorPath = filePrefix + ".connector.circos";
        Files.write(new File(connectorPath).toPath(), createConnectors(data.connectors()));

        String linkPath = filePrefix + ".link.circos";
        Files.write(new File(linkPath).toPath(), createLinks(links));

        String scatterPath = filePrefix + ".scatter.circos";
        Files.write(new File(scatterPath).toPath(), createScatter(segments, links));

        String scatterSglPath = filePrefix + ".scatter.sgl.circos";
        Files.write(new File(scatterSglPath).toPath(), createSglScatter(links));

        String cnaPath = filePrefix + ".cna.circos";
        Files.write(new File(cnaPath).toPath(), createCNA(alterations));

        String mapPath = filePrefix + ".map.circos";
        Files.write(new File(mapPath).toPath(), createMinorAlleleCopyNumber(alterations));

        String centromeresPath = filePrefix + ".centromeres.circos";
        Files.write(new File(centromeresPath).toPath(), highlights(centromereSites));

        String fragileSitesPath = filePrefix + ".fragile.circos";
        Files.write(new File(fragileSitesPath).toPath(), highlights(fragileSites));

        String linePath = filePrefix + ".line_element.circos";
        Files.write(new File(linePath).toPath(), highlights(lineElements));

        String distancesPath = filePrefix + ".distance.circos";
        Files.write(new File(distancesPath).toPath(), createDistances(data.UnadjustedCopyNumbers, alterations));

        String chromosomeRangesPath = filePrefix + ".chromosome.circos";
        Files.write(new File(chromosomeRangesPath).toPath(), createChromosomeRanges(data.ChromosomeRanges));

        String amberBAFsPath = filePrefix + ".amber.circos";
        Files.write(new File(amberBAFsPath).toPath(), createAmberBAFs(data.AmberBAFs));

        String amberPCFPath = filePrefix + ".amber_pcf.circos";
        Files.write(new File(amberPCFPath).toPath(), createAmberPcf(data.PurpleSegments));

        String cobaltRatiosPath = filePrefix + ".cobalt.circos";
        Files.write(new File(cobaltRatiosPath).toPath(), createCobaltRatios(data.CobaltRatios));

        String cobaltPCFPath = filePrefix + ".cobalt_pcf.circos";
        Files.write(new File(cobaltPCFPath).toPath(), createCobaltPcf(data.PurpleSegments));

        return this;
    }

    private List<String> createChromosomeRanges(final List<GenomeRegion> regions)
    {
        List<String> result = Lists.newArrayList();

        for(GenomeRegion region : regions)
        {
            String bandString = new StringJoiner(DELIMITER).add(region.chromosome())
                    .add(String.valueOf(region.start()))
                    .add(String.valueOf(region.end()))
                    .add(ColorPicker.hexContigColor(region.chromosome()))
                    .toString();
            result.add(bandString);
        }

        return result;
    }

    @NotNull
    private List<String> genes(final Map<String, String> geneColours, final List<Gene> genes)
    {
        List<String> result = Lists.newArrayList();
        if(!data.displayGenes())
        {
            return result;
        }

        for(Gene gene : genes)
        {
            String exonString = new StringJoiner(DELIMITER).add(circosContig(gene.chromosome()))
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
    private List<String> exonRank(int totalContigLength, final List<VisGeneExon> exons)
    {
        List<String> result = Lists.newArrayList();
        if(!data.displayGenes())
        {
            return result;
        }

        if(totalContigLength <= MAX_CONTIG_LENGTH_TO_DISPLAY_EXON_RANK)
        {
            for(VisGeneExon exon : exons)
            {
                int position = exon.start() + (exon.end() - exon.start()) / 2;

                String exonString = new StringJoiner(DELIMITER).add(circosContig(exon.chromosome()))
                        .add(String.valueOf(position))
                        .add(String.valueOf(position))
                        .add(String.valueOf(exon.ExonRank))
                        .toString();
                result.add(exonString);

            }
        }

        return result;
    }

    @NotNull
    private List<String> geneName(final List<Gene> genes)
    {
        List<String> result = Lists.newArrayList();
        if(!data.displayGenes())
        {
            return result;
        }

        for(Gene gene : genes)
        {
            StringJoiner exonStringJoiner = new StringJoiner(DELIMITER).add(circosContig(gene.chromosome()))
                    .add(String.valueOf(gene.namePosition()))
                    .add(String.valueOf(gene.namePosition()))
                    .add(gene.name())
                    .add("label_size=" + data.geneLabelSize() + "p");

            result.add(exonStringJoiner.toString());
        }

        return result;
    }

    private List<String> exons(final Map<String, String> geneColours, final List<GenomeRegion> disruptedRegions, final List<VisGeneExon> exons)
    {
        List<String> result = Lists.newArrayList();
        if(!data.displayGenes())
        {
            return result;
        }

        for(VisGeneExon exon : exons)
        {
            String exonString = new StringJoiner(DELIMITER).add(circosContig(exon.chromosome()))
                    .add(String.valueOf(exon.start()))
                    .add(String.valueOf(exon.end()))
                    .add(String.valueOf(1))
                    .add("fill_color=" + geneColours.get(exon.Gene))
                    .toString();
            result.add(exonString);
        }

        for(GenomeRegion disruptedRegion : disruptedRegions)
        {
            String exonString = new StringJoiner(DELIMITER).add(circosContig(disruptedRegion.chromosome()))
                    .add(String.valueOf(disruptedRegion.start()))
                    .add(String.valueOf(disruptedRegion.end()))
                    .add(String.valueOf(1))
                    .add("fill_color=(255,255,255,0.6)")
                    .toString();
            result.add(exonString);
        }

        return result;
    }

    private List<String> createAmberBAFs(List<AmberBAF> amberBAFs)
    {
        List<String> result = Lists.newArrayList();

        for(AmberBAF amberBAF : amberBAFs)
        {
            String amberBAFString =  new StringJoiner(DELIMITER).add(circosContig(amberBAF.Chromosome))
                    .add(String.valueOf(amberBAF.Position))
                    .add(String.valueOf(amberBAF.Position))
                    .add(String.valueOf(amberBAF.TumorBAF))
                    .toString();

            result.add(amberBAFString);
        }

        return result;
    }

    private List<String> createAmberPcf(List<PurpleSegment> purpleSegments)
    {
        List<String> result = Lists.newArrayList();

        for(PurpleSegment purpleSegment : purpleSegments)
        {
            // Ensure histogram mirror images do not cross
            double observedBAF = Math.min(1-purpleSegment.ObservedBAF, purpleSegment.ObservedBAF);

            String purpleSegmentString = new StringJoiner(DELIMITER).add(circosContig(purpleSegment.Chromosome))
                    .add(String.valueOf(purpleSegment.PosStart))
                    .add(String.valueOf(purpleSegment.PosEnd))
                    .add(String.valueOf(observedBAF))
                    .toString();

            result.add(purpleSegmentString);
        }

        return result;
    }

    private List<String> createCobaltRatios(List<CobaltRatio> cobaltRatios)
    {
        List<String> result = Lists.newArrayList();

        for(CobaltRatio cobaltRatio : cobaltRatios)
        {
            String cobaltRatioString =  new StringJoiner(DELIMITER).add(circosContig(cobaltRatio.chromosome()))
                    .add(String.valueOf(cobaltRatio.position()))
                    .add(String.valueOf(cobaltRatio.position()))
                    .add(String.valueOf(cobaltRatio.tumorGCRatio()))
                    .toString();

            result.add(cobaltRatioString);
        }

        return result;
    }

    private List<String> createCobaltPcf(List<PurpleSegment> purpleSegments)
    {
        List<String> result = Lists.newArrayList();

        for(PurpleSegment purpleSegment : purpleSegments)
        {
            String purpleSegmentString = new StringJoiner(DELIMITER).add(circosContig(purpleSegment.Chromosome))
                    .add(String.valueOf(purpleSegment.PosStart))
                    .add(String.valueOf(purpleSegment.PosEnd))
                    .add(String.valueOf(purpleSegment.ObservedTumorRatio))
                    .toString();

            result.add(purpleSegmentString);
        }

        return result;
    }

    private List<String> highlights(final List<GenomeRegion> regions)
    {
        return regions.stream()
                .map(x -> new StringJoiner(DELIMITER).add(circosContig(x.chromosome()))
                        .add(String.valueOf(x.start()))
                        .add(String.valueOf(x.end()))
                        .toString())
                .collect(toList());
    }

    private List<String> createDistances(final List<VisCopyNumber> unadjustedSegment, final List<VisCopyNumber> segments)
    {
        List<String> result = Lists.newArrayList();
        int unadjustedSegments = (int)segments.stream().filter(x -> !x.Truncated).count();
        if(unadjustedSegments <= circosConfig.MaxNumberOfDistanceLabels)
        {
            for(int i = 0; i < unadjustedSegment.size(); i++)
            {
                VisCopyNumber adjusted = segments.get(i);
                VisCopyNumber unadjusted = unadjustedSegment.get(i);
                if(!adjusted.Truncated)
                {
                    String distance = new StringJoiner(DELIMITER).add(circosContig(adjusted.chromosome()))
                            .add(String.valueOf(adjusted.start()))
                            .add(String.valueOf(adjusted.end()))
                            .add(shorthand(unadjusted.bases()))
                            .toString();
                    result.add(distance);
                }
            }

        }

        return result;
    }

    @NotNull
    private List<String> createCNA(final List<VisCopyNumber> alterations)
    {
        List<String> result = Lists.newArrayList();
        for(VisCopyNumber alteration : alterations)
        {
            String cna = new StringJoiner(DELIMITER).add(circosContig(alteration.chromosome()))
                    .add(String.valueOf(alteration.start()))
                    .add(String.valueOf(alteration.end()))
                    .add(String.valueOf(alteration.CopyNumber - 2))
                    .toString();
            result.add(cna);
        }

        return result;
    }

    private List<String> createMinorAlleleCopyNumber(final List<VisCopyNumber> alterations)
    {
        List<String> result = Lists.newArrayList();
        for(VisCopyNumber alteration : alterations)
        {
            String cna = new StringJoiner(DELIMITER).add(circosContig(alteration.chromosome()))
                    .add(String.valueOf(alteration.start()))
                    .add(String.valueOf(alteration.end()))
                    .add(String.valueOf(alteration.minorAlleleCopyNumber() - 1))
                    .toString();
            result.add(cna);
        }

        return result;
    }

    private List<String> createScatter(final List<VisSegment> segments, final List<VisSvData> links)
    {
        int glyphSize = circosConfig.GlyphSize;
        int glyphSizeInner = (int) Math.floor(circosConfig.GlyphSize * 14d / 20d);

        List<String> result = Lists.newArrayList();
        for(VisSegment segment : segments)
        {

            if(segment.Track == 0)
            {
                continue;
            }

            String colorOption = colorPicker.transparentColor(segment.ClusterId, segment.ChainId);
            String startGlyph = scatterGlyph(true, segment, links);
            result.add(scatterEntry(true, segment, colorOption, startGlyph, glyphSize, segment.Frame));
            if(segment.startTerminal() == SegmentTerminal.CENTROMERE)
            {
                result.add(scatterEntry(true, segment, "color=white", startGlyph, glyphSizeInner, 0));
            }

            String endGlyph = scatterGlyph(false, segment, links);
            result.add(scatterEntry(false, segment, colorOption, endGlyph, glyphSize, segment.Frame));
            if(segment.endTerminal() == SegmentTerminal.CENTROMERE)
            {
                result.add(scatterEntry(false, segment, "color=white", endGlyph, glyphSizeInner, 0));
            }
        }

        return result;
    }

    private List<String> createSglScatter(final List<VisSvData> links)
    {
        int glyphSize = circosConfig.GlyphSize;
        int glyphSizeInner = (int) Math.floor(circosConfig.GlyphSize * 14d / 20d);

        List<String> result = Lists.newArrayList();

        // Draw open circles at SGL ends
        for(VisSvData link : links)
        {
            if(link.isValidStart() && !link.isValidEnd())
            {
                String colorOption = colorPicker.transparentColor(link.ClusterId, link.ChainId);
                result.add(scatterSGLEntry(link, colorOption, glyphSize));
                result.add(scatterSGLEntry(link, "color=white", glyphSizeInner));
            }
        }

        return result;
    }

    private String scatterGlyph(boolean isStart, final VisSegment segment, final List<VisSvData> links)
    {
        int location = isStart ? segment.start() : segment.end();
        SegmentTerminal terminal = isStart ? segment.startTerminal() : segment.endTerminal();
        if(terminal != SegmentTerminal.NONE)
        {
            return "square";
        }

        GenomePosition startPosition = GenomePositions.create(segment.chromosome(), location);

        boolean isFoldback =
                VisLinks.findStartLink(startPosition, links).filter(x -> x.InfoStart.equals("FOLDBACK")).isPresent() || VisLinks.findEndLink(
                        startPosition,
                        links).filter(x -> x.InfoEnd.equals("FOLDBACK")).isPresent();

        return isFoldback ? "triangle" : "circle";
    }

    private String scatterEntry(
            boolean isStart, final VisSegment segment, final String color, final String glyph, int glyphSize, int frame)
    {

        int location = isStart ? segment.start() : segment.end();

        return new StringJoiner(DELIMITER).add(circosContig(segment.chromosome()))
                .add(String.valueOf(location))
                .add(String.valueOf(location))
                .add(String.valueOf(segment.Track))
                .add(color + "," + "glyph=" + glyph + ",glyph_size=" + glyphSize + ",frame=" + frame)
                .toString();
    }

    private String scatterSGLEntry(final VisSvData link, final String color, int glyphSize)
    {
        return new StringJoiner(DELIMITER).add(circosContig(link.ChrStart))
                .add(String.valueOf(link.PosStart))
                .add(String.valueOf(link.PosStart))
                .add(String.valueOf(0))
                .add(color + "," + "glyph=circle,glyph_size=" + glyphSize)
                .toString();
    }

    private List<String> createLinks(final List<VisSvData> links)
    {
        final List<String> result = Lists.newArrayList();
        for(final VisSvData link : links)
        {
            if(!HumanChromosome.contains(link.ChrStart) || !HumanChromosome.contains(link.ChrEnd))
                continue;

            final String linkString = new StringJoiner(DELIMITER).add(circosContig(link.ChrStart))
                    .add(String.valueOf(link.PosStart))
                    .add(String.valueOf(link.PosStart))
                    .add(circosContig(link.ChrEnd))
                    .add(String.valueOf(link.PosEnd))
                    .add(String.valueOf(link.PosEnd))
                    .add(colorPicker.transparentColor(link.ClusterId, link.ChainId) + "," + thicknessString(link.JCN) + ",frame="
                            + link.Frame)
                    .toString();
            result.add(linkString);
        }

        return result;
    }

    private List<String> createConnectors(final List<Connector> connectors)
    {
        final List<String> result = Lists.newArrayList();
        for(Connector connector : connectors)
        {
            final double r1 = configWriter.svTrackRelative(connector.track());
            final String start = new StringJoiner(DELIMITER).add(circosContig(connector.chromosome()))
                    .add(String.valueOf(connector.position()))
                    .add(String.valueOf(connector.position()))
                    .add("r1=" + RATIO_FORMAT.format(r1) + "r," + colorPicker.transparentColor(connector.clusterId(), connector.chainId())
                            + "," + thicknessString(connector.ploidy())
                            + ",frame=" + connector.frame()
                    )
                    .toString();
            result.add(start);
        }

        return result;
    }

    @NotNull
    private List<String> createKaryotypes(final Set<GenomePosition> contigLengths)
    {
        final List<String> result = Lists.newArrayList();
        final List<GenomePosition> positions = Lists.newArrayList(contigLengths);
        Collections.sort(positions);

        for(GenomePosition contig : positions)
        {

            final String start = new StringJoiner(" ").add("chr -")
                    .add(circosContig(contig.chromosome()))
                    .add(HumanChromosome.fromString(contig.chromosome()).toString())
                    .add(String.valueOf(1))
                    .add(String.valueOf(Math.max(MIN_KARYOTYPE_LENGTH, contig.position())))
                    .add("chr" + HumanChromosome.fromString(contig.chromosome()).toString())
                    .toString();
            result.add(start);
        }

        return result;
    }

    private List<String> createHistogramTrack(final List<VisSegment> segments)
    {
        final List<String> result = Lists.newArrayList();
        for(final VisSegment segment : segments)
        {
            if(segment.Track > 0)
            {
                double thickness = thicknessPixels(segment.LinkPloidy);
                double r0 = configWriter.svTrackRelative(segment.Track);
                String r0String = "r0=" + RATIO_FORMAT.format(r0) + "r-" + thickness / 2d + "p";
                String r1String = "r1=" + RATIO_FORMAT.format(r0) + "r+" + thickness / 2d + "p";

                final String entry = new StringJoiner(DELIMITER).add(circosContig(segment.chromosome()))
                        .add(String.valueOf(segment.start()))
                        .add(String.valueOf(segment.end()))
                        .add(String.valueOf(segment.Track))
                        .add("fill_" + colorPicker.transparentColor(segment.ClusterId, segment.ChainId) + "," + r0String + ","
                                + r1String + ",frame=" + segment.Frame)
                        .toString();
                result.add(entry);
            }
        }

        return result;
    }

    private List<String> createPositionText(final List<VisSvData> originalLinks, final List<VisSvData> scaledLinks)
    {
        final List<AdjustedPosition> positions = AdjustedPositions.create(originalLinks, scaledLinks);
        if(circosConfig.ExactPosition)
        {
            return createPositionText(1, positions, POSITION_FORMAT::format);
        }

        final List<String> positionsEvery100k = createPositionText(100_000, positions, CircosDataWriter::shorthand);
        if(positionsEvery100k.size() < circosConfig.MaxNumberOfPositionLabels)
        {
            return positionsEvery100k;
        }

        final List<String> positionsEvery1M = createPositionText(1_000_000, positions, CircosDataWriter::shorthand);
        if(positionsEvery1M.size() < circosConfig.MaxNumberOfPositionLabels)
        {
            return positionsEvery1M;
        }

        return createPositionText(10_000_000, positions, CircosDataWriter::shorthand);
    }

    private List<String> createPositionText(int minDistance, final List<AdjustedPosition> positions, final Function<Integer,String> formatter)
    {
        final Set<String> result = Sets.newHashSet();
        final Set<String> contigs = positions.stream().map(GenomePosition::chromosome).collect(Collectors.toSet());

        for(final String contig : contigs)
        {
            int currentPosition = 0;
            for(final AdjustedPosition adjustedPosition : positions)
            {
                if(adjustedPosition.chromosome().equals(contig))
                {
                    int newPosition = adjustedPosition.unadjustedPosition();

                    if(newPosition - minDistance >= currentPosition)
                    {
                        String positionLabel = formatter.apply(adjustedPosition.unadjustedPosition());

                        if(circosConfig.ShowSvId)
                        {
                            positionLabel += String.format(":%d", adjustedPosition.svId());
                        }

                        final String start = new StringJoiner(DELIMITER).add(circosContig(contig))
                                .add(String.valueOf(adjustedPosition.position()))
                                .add(String.valueOf(adjustedPosition.position()))
                                .add(positionLabel)
                                .toString();

                        result.add(start);
                        currentPosition = newPosition;
                    }
                }
            }
        }

        return result.stream().sorted().distinct().collect(toList());
    }

    @NotNull
    private static String circosContig(final String chromosome)
    {
        return "hs" + HumanChromosome.fromString(chromosome);
    }

    @NotNull
    private String thicknessString(double usage)
    {
        return "thickness=" + thicknessPixels(usage);
    }

    private double thicknessPixels(double ploidy)
    {
        return thickness.thicknessPixels(ploidy);
    }

    @NotNull
    static String shorthand(int value)
    {
        if(value < 100)
        {
            return String.valueOf(value);
        }

        if(value < 99_950)
        {
            return String.format(Locale.ENGLISH, "%.1fk", value / 1_000d);
        }

        return String.format(Locale.ENGLISH, "%.1fm", value / 1_000_000d);
    }

}
