package com.hartwig.hmftools.common.version;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.util.List;

import org.apache.logging.log4j.core.util.IOUtils;
import org.jetbrains.annotations.NotNull;

public class VersionInfo {

    private final String resource;

    public VersionInfo(final String resource) {
        this.resource = resource;
    }

    public String version() throws IOException {
        return value("version=", "UNKNOWN");
    }

    public String value(@NotNull final String key, @NotNull final String defaultValue) {
        try {
            for (String entry : readResource().split("\n")) {
                if (entry.startsWith(key)) {
                    return entry.substring(key.length());
                }
            }
        } catch (IOException ignored) {
        }
        return defaultValue;
    }

    public void write(@NotNull final String outputDirectory) throws IOException {
        Files.write(new File(outputDirectory + File.separator + resource).toPath(), read());
    }

    private List<String> read() throws IOException {
        final File file = new File(VersionInfo.class.getClassLoader().getResource(resource).getFile());
        return Files.readAllLines(file.toPath());
    }

    @NotNull
    private String readResource() throws IOException {
        InputStream in = VersionInfo.class.getClassLoader().getResourceAsStream(resource);
        BufferedReader reader = new BufferedReader(new InputStreamReader(in));
        return IOUtils.toString(reader);
    }
}
