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

    private final String sample;
    private final int maxTracks;
    private final String configPath;

    public CircosConfigWriter(@NotNull final String sample, @NotNull final String outputDir, int maxTracks) {
        this.sample = sample;
        this.configPath = outputDir + File.separator + sample + ".circos.conf";
        this.maxTracks = maxTracks;
    }

    public String configPath() {
        return configPath;
    }

    public void writeConfig(int maxTracks, int maxCNATracks) throws IOException {
        final Charset charset = StandardCharsets.UTF_8;
        final String template =
                readResource("/visualisation/cluster.template").replaceAll("SUBSTITUTE_HISTOGRAM", histogramPlots(maxTracks))
                        .replaceAll("SUBSTITUTE_MAX_CNA_INV", String.valueOf(1d / maxCNATracks))
                        .replaceAll("SUBSTITUTE_MAX_CNA", String.valueOf(maxCNATracks))
                        .replaceAll("SUBSTITUTE_MAX_INV", String.valueOf(1d / maxTracks))
                        .replaceAll("SUBSTITUTE_MAX", String.valueOf(maxTracks))
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
