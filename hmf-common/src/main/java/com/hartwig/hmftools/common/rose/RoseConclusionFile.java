package com.hartwig.hmftools.common.rose;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class RoseConclusionFile {
    private static final Logger LOGGER = LogManager.getLogger(RoseConclusionFile.class);

    private static final String EXTENSION = ".rose.tsv";

    private RoseConclusionFile() {
    }

    @NotNull
    public static String generateFilename(@NotNull String basePath, @NotNull String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    @NotNull
    public static String read(@NotNull String file) throws IOException {
        LOGGER.info("Reading ROSE clinical conslusion {}", file);
        String conclusion = Strings.EMPTY;
        List<String> lines = Files.readAllLines(new File(file).toPath());

        for (String line : lines) {
            conclusion = conclusion.concat(line + " \n ");
        }
        return conclusion;
    }


    public static void write(@NotNull String file, @NotNull ActionabilityConclusion actionabilityConclusion)
            throws IOException {

        List<String> lines = Lists.newArrayList();
        for (String conclusion : actionabilityConclusion.conclusion()) {
            lines.add(conclusion);
        }
        Files.write(new File(file).toPath(), lines);
    }
}