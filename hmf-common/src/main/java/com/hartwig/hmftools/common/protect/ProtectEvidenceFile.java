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
                .add("transcript")
                .add("isCanonical")
                .add("event")
                .add("evidenceType")
                .add("eventIsHighDriver")
                .add("rangeRank")
                .add("germline")
                .add("reported")
                .add("treatment")
                .add("onLabel")
                .add("level")
                .add("direction")
                .add("evidenceUrls")
                .add("sources")
                .add("sourceEvent")
                .add("sourceUrls")
                .toString();
    }

    @NotNull
    private static String toLine(@NotNull ProtectEvidence evidence) {
        StringJoiner evidenceUrlJoiner = new StringJoiner(SUBFIELD_DELIMITER);
        for (String url : evidence.evidenceUrls()) {
            evidenceUrlJoiner.add(url);
        }

        StringJoiner sourceUrlJoiner = new StringJoiner(SUBFIELD_DELIMITER);
        for (String url : evidence.sourceUrls()) {
            sourceUrlJoiner.add(url);
        }

        StringJoiner sourceJoiner = new StringJoiner(SUBFIELD_DELIMITER);
        for (Knowledgebase source : evidence.sources()) {
            sourceJoiner.add(source.technicalDisplay());
        }

        return new StringJoiner(FIELD_DELIMITER).add(nullToEmpty(evidence.gene()))
                .add(evidence.transcript())
                .add(String.valueOf(evidence.isCanonical()))
                .add(evidence.event())
                .add(evidence.evidenceType().toString())
                .add(nullToEmpty(evidence.eventIsHighDriver()))
                .add(nullToEmpty(evidence.rangeRank()))
                .add(String.valueOf(evidence.germline()))
                .add(String.valueOf(evidence.reported()))
                .add(evidence.treatment())
                .add(String.valueOf(evidence.onLabel()))
                .add(evidence.level().toString())
                .add(evidence.direction().toString())
                .add(evidenceUrlJoiner.toString())
                .add(sourceJoiner.toString())
                .add(evidence.sourceEvent())
                .add(sourceUrlJoiner.toString())
                .toString();
    }

    @NotNull
    private static String nullToEmpty(@Nullable Boolean booleanValue) {
        return booleanValue != null ? Boolean.toString(booleanValue) : Strings.EMPTY;
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

        String evidenceUrlField = values[fields.get("evidenceUrls")];
        Set<String> evidenceUrlurls = !evidenceUrlField.isEmpty() ? Sets.newHashSet(evidenceUrlField.split(SUBFIELD_DELIMITER)) : Sets.newHashSet();

        String sourceEventField = values[fields.get("sourceEvent")];
        String sourceEvent = !sourceEventField.isEmpty() ? sourceEventField : Strings.EMPTY;

        String eventIsHighDriverField = values[fields.get("eventIsHighDriver")];
        Boolean eventIsHighDriver = !eventIsHighDriverField.isEmpty() ? Boolean.parseBoolean(eventIsHighDriverField) : null;

        String sourceUrlField = values[fields.get("sourceUrls")];
        Set<String> sourceUrlurls = !sourceUrlField.isEmpty() ? Sets.newHashSet(sourceUrlField.split(SUBFIELD_DELIMITER)) : Sets.newHashSet();

        String transcriptField = values[fields.get("transcript")];
        String transcript = !transcriptField.isEmpty() ? transcriptField : Strings.EMPTY;

        return ImmutableProtectEvidence.builder()
                .gene(emptyToNullString(values[fields.get("gene")]))
                .transcript(transcript)
                .isCanonical(Boolean.parseBoolean(values[fields.get("isCanonical")]))
                .event(values[fields.get("event")])
                .eventIsHighDriver(eventIsHighDriver)
                .evidenceType(ProtectEvidenceType.valueOf(values[fields.get("evidenceType")]))
                .rangeRank(emptyToNullInteger(values[fields.get("rangeRank")]))
                .germline(Boolean.parseBoolean(values[fields.get("germline")]))
                .reported(Boolean.parseBoolean(values[fields.get("reported")]))
                .treatment(values[fields.get("treatment")])
                .onLabel(Boolean.parseBoolean(values[fields.get("onLabel")]))
                .level(EvidenceLevel.valueOf(values[fields.get("level")]))
                .direction(EvidenceDirection.valueOf(values[fields.get("direction")]))
                .evidenceUrls(evidenceUrlurls)
                .sources(Knowledgebase.fromCommaSeparatedTechnicalDisplayString(values[fields.get("sources")]))
                .sourceEvent(sourceEvent)
                .sourceUrls(sourceUrlurls)
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
