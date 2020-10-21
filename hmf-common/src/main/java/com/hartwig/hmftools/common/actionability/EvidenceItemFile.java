package com.hartwig.hmftools.common.actionability;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class EvidenceItemFile {

    private static final String DELIMITER = "\t";
    private static final String EXTENSION = ".evidence.tsv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static void write(@NotNull String file, @NotNull List<EvidenceItem> evidence) throws IOException {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        lines.addAll(evidence.stream().map(EvidenceItemFile::toLine).collect(Collectors.toList()));
        Files.write(new File(file).toPath(), lines);
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER).add("event")
                .add("source")
                .add("reference")
                .add("drug")
                .add("drugsType")
                .add("level")
                .add("response")
                .add("onLabel")
                .add("cancerType")
                .add("scope")
                .toString();
    }

    @NotNull
    private static String toLine(@NotNull EvidenceItem evidence) {
        return new StringJoiner(DELIMITER).add(evidence.event())
                .add(evidence.source().toString())
                .add(evidence.reference())
                .add(evidence.drug())
                .add(evidence.drugsType())
                .add(evidence.level().toString())
                .add(evidence.response())
                .add(String.valueOf(evidence.isOnLabel()))
                .add(evidence.cancerType())
                .add(evidence.scope().toString())
                .toString();
    }

}
