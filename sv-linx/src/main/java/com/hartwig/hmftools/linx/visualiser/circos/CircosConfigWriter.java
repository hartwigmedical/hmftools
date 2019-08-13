package com.hartwig.hmftools.linx.visualiser.circos;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.StringJoiner;
import java.util.function.Function;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.linx.visualiser.SvCircosConfig;

import org.apache.logging.log4j.core.util.IOUtils;
import org.jetbrains.annotations.NotNull;

public class CircosConfigWriter
{

    static final double PIXELS = 1500;

    //    private static final double EXON_RANK_INNER_RADIUS = 0.95;

    //    private static final double EXON_OUTER_RADIUS = 0.95;
    static final double EXON_INNER_RADIUS = 0.9;

    private static final double GENE_OUTER_RADIUS = 0.95;
    private static final double GENE_INNER_RADIUS = 0.93;
    //
    //    private static final double SEGMENT_OUTER_RADIUS = 0.875;
    //    private static final double SEGMENT_INNER_RADIUS = 0.5;
    //
    //    private static final double MAP_OUTER_RADIUS = 0.275;
    //        private static final double INNER_RADIUS = 0.175;
    //
    //    private static final double CNA_OUTER_RADIUS = 0.475;
    //    private static final double CNA_INNER_RADIUS = 0.3;

    private final String sample;
    private final String configPath;
    private final String outputDir;
    private final CircosData circosData;
    private final SvCircosConfig config;

    private final double exonOuterRadius;
    private final double exonInnerRadius;

    private final double segmentOuterRadius;
    private final double segmentInnerRadius;

    private final double copyNumberOuterRadius;
    private final double copyNumberMiddleRadius;
    private final double copyNumberInnerRadius;

    private final double mapOuterRadius;
    private final double mapMiddleRadius;
    private final double mapInnerRadius;

    public CircosConfigWriter(@NotNull final String sample, @NotNull final String outputDir, @NotNull final CircosData data,
            @NotNull final SvCircosConfig config)
    {
        this.sample = sample;
        this.configPath = outputDir + File.separator + sample + ".circos.conf";
        this.circosData = data;
        this.config = config;
        this.outputDir = outputDir;

        double gapSize = config.gapSize();
        double geneRelativeSize = data.exons().isEmpty() ? 0 : config.geneRelativeSize();
        double segmentRelativeSize = config.segmentRelativeSize();
        double copyNumberRelativeSize = config.copyNumberRelativeSize();

        double totalRelativeSize = geneRelativeSize + segmentRelativeSize + copyNumberRelativeSize;

        boolean displayGenes = Doubles.greaterThan(geneRelativeSize, 0);
        int numberOfGaps = displayGenes ? 4 : 3;

        double totalSpaceAvailable = 1 - numberOfGaps * gapSize - config.innerRadius();
        double purpleSpaceAvailable = copyNumberRelativeSize / totalRelativeSize * totalSpaceAvailable;
        int cnaGainTracks = Math.max(2, (int) Math.round(Math.ceil(data.maxCopyNumber() - 2)));
        int mapGainTracks = Math.max(1, (int) Math.round(Math.ceil(data.maxMinorAllelePloidy() - 1)));
        double purpleTrackSize = purpleSpaceAvailable / (1 + 2 + cnaGainTracks + mapGainTracks);

        if (displayGenes)
        {
            exonOuterRadius = 1 - gapSize;
            exonInnerRadius = exonOuterRadius - geneRelativeSize / totalRelativeSize * totalSpaceAvailable;
            segmentOuterRadius = exonInnerRadius - gapSize;
        }
        else
        {
            exonOuterRadius = 0;
            exonInnerRadius = 0;
            segmentOuterRadius = 1 - gapSize;
        }
        segmentInnerRadius = segmentOuterRadius - segmentRelativeSize / totalRelativeSize * totalSpaceAvailable;

        copyNumberOuterRadius = segmentInnerRadius - gapSize;
        copyNumberMiddleRadius = copyNumberOuterRadius - cnaGainTracks * purpleTrackSize;
        copyNumberInnerRadius = copyNumberMiddleRadius - 2 * purpleTrackSize;

        mapOuterRadius = copyNumberInnerRadius - gapSize;
        mapMiddleRadius = mapOuterRadius - mapGainTracks * purpleTrackSize;
        mapInnerRadius = mapMiddleRadius - 1 * purpleTrackSize;
    }

    public String configPath()
    {
        return configPath;
    }

    public double svTrackRelative(int track)
    {
        double difference = segmentOuterRadius() - segmentInnerRadius;
        double singleTrack = difference / circosData.maxTracks();

        return segmentInnerRadius + track * singleTrack;
    }

    private double segmentOuterRadius()
    {
        return segmentOuterRadius;
    }

    public void writeConfig()
            throws IOException
    {
        int chromosomeCount = circosData.contigLengths().size();
        int totalContigLength = circosData.totalContigLength();

        int cnaMaxTracks = Math.max(2, (int) Math.round(Math.ceil(circosData.maxCopyNumber() - 2)));

        int mapMaxTracks = Math.max(1, (int) Math.round(Math.ceil(circosData.maxMinorAllelePloidy() - 1)));
        final Charset charset = StandardCharsets.UTF_8;
        final String template =
                readResource("/visualisation/cluster.template")
                        .replaceAll("SUBSTITUTE_IDEOGRAM_RADIUS", config.displayPosition() ? "0.90" : "0.99")
                        .replaceAll("SUBSTITUTE_IDEOGRAM_SPACING", chromosomeCount > 1 ? "0.005r" : 0.005 * totalContigLength + "u")

                        .replaceAll("SUBSTITUTE_EXON_RANK_INNER_RADIUS", String.valueOf(exonInnerRadius))

                        .replaceAll("SUBSTITUTE_EXON_INNER_RADIUS", String.valueOf(exonInnerRadius))
                        .replaceAll("SUBSTITUTE_EXON_OUTER_RADIUS", String.valueOf(exonOuterRadius))
                        .replaceAll("SUBSTITUTE_GENE_INNER_RADIUS", String.valueOf(GENE_INNER_RADIUS))
                        .replaceAll("SUBSTITUTE_GENE_OUTER_RADIUS", String.valueOf(GENE_OUTER_RADIUS))

                        .replaceAll("SUBSTITUTE_SV_INNER_RADIUS", String.valueOf(segmentInnerRadius))
                        .replaceAll("SUBSTITUTE_SV_OUTER_RADIUS", String.valueOf(segmentOuterRadius))

                        .replaceAll("SUBSTITUTE_MAP_INNER_RADIUS", String.valueOf(mapInnerRadius))
                        .replaceAll("SUBSTITUTE_MAP_OUTER_RADIUS", String.valueOf(mapOuterRadius))
                        .replaceAll("SUBSTITUTE_MAP_MIDDLE_RADIUS", String.valueOf(mapMiddleRadius))
                        .replaceAll("SUBSTITUTE_MAP_GAIN_MAX", String.valueOf(mapMaxTracks))
                        .replaceAll("SUBSTITUTE_MAP_GAIN_SPACING", String.valueOf(1d / mapMaxTracks))

                        .replaceAll("SUBSTITUTE_CNA_INNER_RADIUS", String.valueOf(copyNumberInnerRadius))
                        .replaceAll("SUBSTITUTE_CNA_OUTER_RADIUS", String.valueOf(copyNumberOuterRadius))
                        .replaceAll("SUBSTITUTE_CNA_MIDDLE_RADIUS", String.valueOf(copyNumberMiddleRadius))
                        .replaceAll("SUBSTITUTE_CNA_GAIN_MAX", String.valueOf(cnaMaxTracks))
                        .replaceAll("SUBSTITUTE_CNA_GAIN_AXIS_POSITION", cnaAxisPositions(cnaMaxTracks))

                        .replaceAll("SUBSTITUTE_SV_SPACING", String.valueOf(1d / circosData.maxTracks()))
                        .replaceAll("SUBSTITUTE_SV_MAX", String.valueOf(circosData.maxTracks()))

                        .replaceAll("SUBSTITUTE_SAMPLE", sample);

        Files.write(new File(configPath).toPath(), template.getBytes(charset));
    }

    public void writeCytobands() throws IOException
    {
        final String template = readResource("/r/cytoBand.txt");
        Files.write(new File(outputDir + File.separator + "cytoBand.txt").toPath(), template.getBytes(StandardCharsets.UTF_8));
    }

    @NotNull
    private String readResource(@NotNull final String resource) throws IOException
    {
        InputStream in = getClass().getResourceAsStream(resource);
        BufferedReader reader = new BufferedReader(new InputStreamReader(in));
        return IOUtils.toString(reader);
    }

    @NotNull
    static String cnaAxisPositions(int maxTracks)
    {
        StringJoiner builder = new StringJoiner(",");

        final double rel = 1d / maxTracks;
        final Function<Integer, String> relString = i -> String.valueOf(Math.round(i * rel * 10000) / 10000d) + "r";

        for (int i = 1; i <= Math.min(7, maxTracks); i++)
        {
            builder.add(relString.apply(i));
        }

        for (int i = 8; i <= maxTracks; i += 10)
        {
            builder.add(relString.apply(i));
        }

        return builder.toString();
    }

}
