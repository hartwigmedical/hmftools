package com.hartwig.hmftools.common.version;

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

public class VersionInfo {

    @NotNull
    private final String resource;

    public VersionInfo(@NotNull final String resource) {
        this.resource = resource;
    }

    @NotNull
    public String version() throws IOException {
        return value("version=", "UNKNOWN");
    }

    @NotNull
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
        final String content = readResource();
        final Charset charset = StandardCharsets.UTF_8;
        Files.write(new File(outputDirectory + File.separator + resource).toPath(), content.getBytes(charset));
    }

    @NotNull
    private String readResource() throws IOException {
        InputStream in = VersionInfo.class.getClassLoader().getResourceAsStream(resource);
        BufferedReader reader = new BufferedReader(new InputStreamReader(in));
        return IOUtils.toString(reader);
    }
}
