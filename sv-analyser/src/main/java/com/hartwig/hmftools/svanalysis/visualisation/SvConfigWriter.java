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

public class SvConfigWriter {

    private final String sample;
    private final String outputDir;

    public SvConfigWriter(final String sample, final String outputDir) {
        this.sample = sample;
        this.outputDir = outputDir;
    }

    public void writeConfig(int levels) throws IOException {
        final Charset charset = StandardCharsets.UTF_8;
        final String template = readResource("/visualisation/cluster.template")
                .replaceAll("SUBSTITUTE_HISTOGRAM", histogramPlots(levels))
                .replaceAll("SUBSTITUTE_MAX", String.valueOf(levels))
                .replaceAll("SUBSTITUTE_SAMPLE", sample);

        Files.write(new File("/Users/jon/hmf/analysis/sv/SvWriter.circos.conf").toPath(), template.getBytes(charset));

    }

    @NotNull
    private String histogramPlots(int levels) throws IOException {
        final StringBuilder builder = new StringBuilder();

        final String histogramTemplate = readResource("/visualisation/cluster.template.histogram");
        for (int i = 1; i <= levels; i++) {
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
