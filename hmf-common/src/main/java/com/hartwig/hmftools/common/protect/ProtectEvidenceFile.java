package com.hartwig.hmftools.common.protect;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.jetbrains.annotations.NotNull;

public final class ProtectEvidenceFile {

    private static final String EXTENSION = ".protect.tsv";
    private static final String FIELD_DELIMITER = "\t";

    private static final String SUBFIELD_DELIMITER = ",";

    private ProtectEvidenceFile() {
    }

    @NotNull
    public static String generateFilename(@NotNull String basePath, @NotNull String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static void write(@NotNull String file, @NotNull List<ProtectEvidence> evidence) throws IOException {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        lines.addAll(evidence.stream().map(ProtectEvidenceFile::toLine).collect(Collectors.toList()));
        Files.write(new File(file).toPath(), lines);
    }

    @NotNull
    public static List<ProtectEvidence> read(@NotNull String file) throws IOException {
        List<ProtectEvidence> evidence = Lists.newArrayList();
        List<String> lines = Files.readAllLines(new File(file).toPath());
        // Skip header
        for (String line : lines.subList(1, lines.size())) {
            evidence.add(fromLine(line));
        }
        return evidence;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(FIELD_DELIMITER).add("event")
                .add("germline")
                .add("reported")
                .add("treatment")
                .add("onLabel")
                .add("level")
                .add("direction")
                .add("sources")
                .add("urls")
                .toString();
    }

    @NotNull
    private static String toLine(@NotNull ProtectEvidence evidence) {
        StringJoiner urlJoiner = new StringJoiner(SUBFIELD_DELIMITER);
        for (String url : evidence.urls()) {
            urlJoiner.add(url);
        }

        StringJoiner sourceJoiner = new StringJoiner(SUBFIELD_DELIMITER);
        for (Knowledgebase source : evidence.sources()) {
            sourceJoiner.add(source.technicalDisplay());
        }

        return new StringJoiner(FIELD_DELIMITER).add(evidence.genomicEvent())
                .add(String.valueOf(evidence.germline()))
                .add(String.valueOf(evidence.reported()))
                .add(evidence.treatment())
                .add(String.valueOf(evidence.onLabel()))
                .add(evidence.level().toString())
                .add(evidence.direction().toString())
                .add(sourceJoiner.toString())
                .add(urlJoiner.toString())
                .toString();
    }

    @NotNull
    private static ProtectEvidence fromLine(@NotNull String line) {
        String[] values = line.split(FIELD_DELIMITER);

        Set<String> urls = values.length == 9 ? Sets.newHashSet(values[8].split(SUBFIELD_DELIMITER)) : Sets.newHashSet();
        return ImmutableProtectEvidence.builder().genomicEvent(values[0])
                .germline(Boolean.parseBoolean(values[1]))
                .reported(Boolean.parseBoolean(values[2]))
                .treatment(values[3])
                .onLabel(Boolean.parseBoolean(values[4]))
                .level(EvidenceLevel.valueOf(values[5]))
                .direction(EvidenceDirection.valueOf(values[6]))
                .sources(Knowledgebase.fromCommaSeparatedTechnicalDisplayString(values[7]))
                .urls(urls)
                .build();
    }
}
