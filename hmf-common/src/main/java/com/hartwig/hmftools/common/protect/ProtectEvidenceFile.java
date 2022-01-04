package com.hartwig.hmftools.common.protect;

import static com.hartwig.hmftools.common.sv.linx.LinxCluster.DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.utils.FileWriterUtils;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

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

        Map<String, Integer> fields = FileWriterUtils.createFieldsIndexMap(lines.get(0), DELIMITER);
        for (String line : lines.subList(1, lines.size())) {
            evidence.add(fromLine(fields, line));
        }
        return evidence;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(FIELD_DELIMITER).add("gene")
                .add("event")
                .add("evidenceType")
                .add("rangeRank")
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

        return new StringJoiner(FIELD_DELIMITER).add(nullToEmpty(evidence.gene()))
                .add(evidence.event())
                .add(evidence.evidenceType().toString())
                .add(nullToEmpty(evidence.rangeRank()))
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
    private static String nullToEmpty(@Nullable String string) {
        return string != null ? string : Strings.EMPTY;
    }

    @NotNull
    private static String nullToEmpty(@Nullable Integer integer) {
        return integer != null ? String.valueOf(integer) : Strings.EMPTY;
    }

    @NotNull
    private static ProtectEvidence fromLine(@NotNull Map<String, Integer> fields, @NotNull String line) {
        String[] values = line.split(FIELD_DELIMITER, -1);

        String urlField = values[fields.get("urls")];
        Set<String> urls = !urlField.isEmpty() ? Sets.newHashSet(urlField.split(SUBFIELD_DELIMITER)) : Sets.newHashSet();

        return ImmutableProtectEvidence.builder()
                .gene(emptyToNullString(values[fields.get("gene")]))
                .event(values[fields.get("event")])
                .evidenceType(ProtectEvidenceType.valueOf(values[fields.get("evidenceType")]))
                .rangeRank(emptyToNullInteger(values[fields.get("rangeRank")]))
                .germline(Boolean.parseBoolean(values[fields.get("germline")]))
                .reported(Boolean.parseBoolean(values[fields.get("reported")]))
                .treatment(values[fields.get("treatment")])
                .onLabel(Boolean.parseBoolean(values[fields.get("onLabel")]))
                .level(EvidenceLevel.valueOf(values[fields.get("level")]))
                .direction(EvidenceDirection.valueOf(values[fields.get("direction")]))
                .sources(Knowledgebase.fromCommaSeparatedTechnicalDisplayString(values[fields.get("sources")]))
                .urls(urls)
                .build();
    }

    @Nullable
    private static String emptyToNullString(@NotNull String value) {
        return !value.isEmpty() ? value : null;
    }

    @Nullable
    private static Integer emptyToNullInteger(@NotNull String value) {
        return !value.isEmpty() ? Integer.valueOf(value) : null;
    }
}
