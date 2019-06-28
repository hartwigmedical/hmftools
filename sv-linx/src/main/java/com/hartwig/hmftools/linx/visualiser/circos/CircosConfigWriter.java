package com.hartwig.hmftools.linx.visualiser.circos;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import org.apache.logging.log4j.core.util.IOUtils;
import org.jetbrains.annotations.NotNull;

public class CircosConfigWriter
{

    static final double PIXELS = 1500;

    private static final double EXON_RANK_INNER_RADIUS = 0.95;

    private static final double EXON_OUTER_RADIUS = 0.95;
    static final double EXON_INNER_RADIUS = 0.9;

    private static final double GENE_OUTER_RADIUS = 0.95;
    private static final double GENE_INNER_RADIUS = 0.93;

    private static final double SEGMENT_OUTER_RADIUS = 0.875;
    private static final double SEGMENT_INNER_RADIUS = 0.5;

    private static final double MAP_OUTER_RADIUS = 0.275;
    private static final double MAP_INNER_RADIUS = 0.175;

    private static final double CNA_OUTER_RADIUS = 0.475;
    private static final double CNA_INNER_RADIUS = 0.3;

    private final String sample;
    private final String configPath;
    private final CircosData circosData;

    public CircosConfigWriter(@NotNull final String sample, @NotNull final String outputDir, @NotNull final CircosData circosData)
    {
        this.sample = sample;
        this.configPath = outputDir + File.separator + sample + ".circos.conf";
        this.circosData = circosData;
    }

    public String configPath()
    {
        return configPath;
    }

    public double svTrackRelative(int track)
    {
        double difference = segmentOuterRadius() - SEGMENT_INNER_RADIUS;
        double singleTrack = difference / circosData.maxTracks();

        return SEGMENT_INNER_RADIUS + track * singleTrack;
    }

    private double segmentOuterRadius()
    {
        return circosData.displayGenes() ? SEGMENT_OUTER_RADIUS : EXON_OUTER_RADIUS;
    }

    public void writeConfig()
            throws IOException
    {
        int chromosomeCount = circosData.contigLengths().size();
        int totalContigLength = circosData.contigLengths().values().stream().mapToInt(x -> x).sum();

        int cnaMaxTracks = Math.max(2, (int) Math.round(Math.ceil(circosData.maxCopyNumber() - 2)));
        double cnaMiddleRadius = CNA_INNER_RADIUS + 2 * (CNA_OUTER_RADIUS - CNA_INNER_RADIUS) / (cnaMaxTracks + 2);

        int mapMaxTracks = Math.max(1, (int) Math.round(Math.ceil(circosData.maxMinorAllelePloidy() - 1)));
        double mapMiddleRadius = MAP_INNER_RADIUS + (MAP_OUTER_RADIUS - MAP_INNER_RADIUS) / (mapMaxTracks + 1);

        final Charset charset = StandardCharsets.UTF_8;
        final String template =
                readResource("/visualisation/cluster.template")
                        .replaceAll("SUBSTITUTE_IDEOGRAM_SPACING", chromosomeCount > 1 ? "0.005r" : 0.005 * totalContigLength + "u")

                        .replaceAll("SUBSTITUTE_EXON_RANK_INNER_RADIUS", String.valueOf(EXON_RANK_INNER_RADIUS))

                        .replaceAll("SUBSTITUTE_EXON_INNER_RADIUS", String.valueOf(EXON_INNER_RADIUS))
                        .replaceAll("SUBSTITUTE_EXON_OUTER_RADIUS", String.valueOf(EXON_OUTER_RADIUS))
                        .replaceAll("SUBSTITUTE_GENE_INNER_RADIUS", String.valueOf(GENE_INNER_RADIUS))
                        .replaceAll("SUBSTITUTE_GENE_OUTER_RADIUS", String.valueOf(GENE_OUTER_RADIUS))

                        .replaceAll("SUBSTITUTE_SV_INNER_RADIUS", String.valueOf(SEGMENT_INNER_RADIUS))
                        .replaceAll("SUBSTITUTE_SV_OUTER_RADIUS", String.valueOf(segmentOuterRadius()))

                        .replaceAll("SUBSTITUTE_MAP_INNER_RADIUS", String.valueOf(MAP_INNER_RADIUS))
                        .replaceAll("SUBSTITUTE_MAP_OUTER_RADIUS", String.valueOf(MAP_OUTER_RADIUS))
                        .replaceAll("SUBSTITUTE_MAP_MIDDLE_RADIUS", String.valueOf(mapMiddleRadius))
                        .replaceAll("SUBSTITUTE_MAP_GAIN_MAX", String.valueOf(mapMaxTracks))
                        .replaceAll("SUBSTITUTE_MAP_GAIN_SPACING", String.valueOf(1d / mapMaxTracks))

                        .replaceAll("SUBSTITUTE_CNA_INNER_RADIUS", String.valueOf(CNA_INNER_RADIUS))
                        .replaceAll("SUBSTITUTE_CNA_OUTER_RADIUS", String.valueOf(CNA_OUTER_RADIUS))
                        .replaceAll("SUBSTITUTE_CNA_MIDDLE_RADIUS", String.valueOf(cnaMiddleRadius))
                        .replaceAll("SUBSTITUTE_CNA_GAIN_MAX", String.valueOf(cnaMaxTracks))
                        .replaceAll("SUBSTITUTE_CNA_GAIN_SPACING", String.valueOf(1d / cnaMaxTracks))

                        .replaceAll("SUBSTITUTE_SV_SPACING", String.valueOf(1d / circosData.maxTracks()))
                        .replaceAll("SUBSTITUTE_SV_MAX", String.valueOf(circosData.maxTracks()))

                        .replaceAll("SUBSTITUTE_SAMPLE", sample);

        Files.write(new File(configPath).toPath(), template.getBytes(charset));
    }

    @NotNull
    private String readResource(@NotNull final String resource) throws IOException
    {
        InputStream in = getClass().getResourceAsStream(resource);
        BufferedReader reader = new BufferedReader(new InputStreamReader(in));
        return IOUtils.toString(reader);
    }
}
