package com.hartwig.hmftools.common.protect;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class ProtectEvidenceItemFile {

    private static final String DELIMITER = "\t";
    private static final String EXTENSION = ".protect.tsv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static void write(@NotNull String file, @NotNull List<ProtectEvidenceItem> evidence) throws IOException {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        lines.addAll(evidence.stream().map(ProtectEvidenceItemFile::toLine).collect(Collectors.toList()));
        Files.write(new File(file).toPath(), lines);
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER).add("event")
                .add("source")
                .add("treatment")
                .add("onLabel")
                .add("level")
                .add("direction")
                .add("reported")
                .add("url")
                .toString();
    }

    @NotNull
    private static String toLine(@NotNull ProtectEvidenceItem evidence) {
        return new StringJoiner(DELIMITER).add(evidence.genomicEvent())
                .add(evidence.source().toString())
                .add(evidence.treatment())
                .add(String.valueOf(evidence.onLabel()))
                .add(evidence.level().toString())
                .add(evidence.direction().toString())
                .add(String.valueOf(evidence.reported()))
                .add(evidence.url())
                .toString();
    }

}
