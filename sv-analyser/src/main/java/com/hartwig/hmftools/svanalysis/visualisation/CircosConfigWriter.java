package com.hartwig.hmftools.svanalysis.visualisation;

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

public class CircosConfigWriter {

    private static final double RADIUS_PIXELS = (0.9 * 1500 - 50);

    private static final double SEGMENT_INNER_RADIUS = 0.5;
    private static final double SEGMENT_OUTER_RADIUS = 0.975;

    private static final double CNA_INNER_RADIUS = 0.3;
    private static final double CNA_OUTER_RADIUS = 0.475;

    private final String sample;
    private final String configPath;

    public CircosConfigWriter(@NotNull final String sample, @NotNull final String outputDir) {
        this.sample = sample;
        this.configPath = outputDir + File.separator + sample + ".circos.conf";
    }

    public String configPath() {
        return configPath;
    }

    public static double svTrackPixels(int maxTracks, int track) {
        double start = SEGMENT_INNER_RADIUS * RADIUS_PIXELS;
        double end = SEGMENT_OUTER_RADIUS * RADIUS_PIXELS;

        double difference = end - start;
        double singleTrack = difference / maxTracks;

        return start + track * singleTrack;
    }

    public void writeConfig(int maxTracks, final double maxCopyNumber) throws IOException {

        int cnaMaxTracks = Math.max(2, (int) Math.round(Math.ceil(maxCopyNumber - 2)));
        double cnaMiddleRadius = CNA_INNER_RADIUS + 2 * (CNA_OUTER_RADIUS - CNA_INNER_RADIUS) / (cnaMaxTracks + 2);

        final Charset charset = StandardCharsets.UTF_8;
        final String template =
                readResource("/visualisation/cluster.template").replaceAll("SUBSTITUTE_HISTOGRAM", histogramPlots(maxTracks))

                        .replaceAll("SUBSTITUTE_SV_INNER_RADIUS", String.valueOf(SEGMENT_INNER_RADIUS))
                        .replaceAll("SUBSTITUTE_SV_OUTER_RADIUS", String.valueOf(SEGMENT_OUTER_RADIUS))

                        .replaceAll("SUBSTITUTE_CNA_INNER_RADIUS", String.valueOf(CNA_INNER_RADIUS))
                        .replaceAll("SUBSTITUTE_CNA_OUTER_RADIUS", String.valueOf(CNA_OUTER_RADIUS))
                        .replaceAll("SUBSTITUTE_CNA_MIDDLE_RADIUS", String.valueOf(cnaMiddleRadius))
                        .replaceAll("SUBSTITUTE_CNA_GAIN_MAX", String.valueOf(cnaMaxTracks))
                        .replaceAll("SUBSTITUTE_CNA_GAIN_SPACING", String.valueOf(1d / cnaMaxTracks))

                        .replaceAll("SUBSTITUTE_SV_SPACING", String.valueOf(1d / maxTracks))
                        .replaceAll("SUBSTITUTE_SV_MAX", String.valueOf(maxTracks))

                        .replaceAll("SUBSTITUTE_SAMPLE", sample);

        Files.write(new File(configPath).toPath(), template.getBytes(charset));
    }

    @NotNull
    private String histogramPlots(int maxTracks) throws IOException {
        final StringBuilder builder = new StringBuilder();

        final String histogramTemplate = readResource("/visualisation/cluster.template.histogram");
        for (int i = 1; i <= maxTracks; i++) {
            final String level = histogramTemplate.replaceAll("SUBSTITUTE_CONDITION", String.valueOf(i));
            builder.append(level);
        }

        return builder.toString();
    }

    @NotNull
    private String readResource(@NotNull final String resource) throws IOException {
        InputStream in = getClass().getResourceAsStream(resource);
        BufferedReader reader = new BufferedReader(new InputStreamReader(in));
        return IOUtils.toString(reader);
    }
}
