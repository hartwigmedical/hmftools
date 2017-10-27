package com.hartwig.hmftools.purple;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import org.jetbrains.annotations.NotNull;

class PurityPloidyEstimateVersion {

    static String version() throws IOException {
        final String versionString = read().get(0);
        return versionString.substring(versionString.indexOf("=") + 1);
    }

    static void write(@NotNull final String outputDirectory) throws IOException {
        Files.write(new File(outputDirectory + File.separator + "purple.version").toPath(), read());
    }

    private static List<String> read() throws IOException {
        final File file = new File(PurityPloidyEstimateVersion.class.getClassLoader().getResource("purple.version").getFile());
        return Files.readAllLines(file.toPath());
    }
}
