package com.hartwig.hmftools.common.rose;

import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class RoseConclusionFile {

    private static final String EXTENSION = ".rose.tsv";
    private static final String FIELD_DELIMITER = "\t";

    private RoseConclusionFile() {
    }

    @NotNull
    public static String generateFilename(@NotNull String basePath, @NotNull String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    @NotNull
    public static String read(@NotNull String file) throws IOException {
        String conclusion = Strings.EMPTY;
        List<String> lines = Files.readAllLines(new File(file).toPath());

        for (String line : lines) {
            conclusion = conclusion.concat(line + " <enter> ");
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