package com.hartwig.hmftools.common.protect;

import static com.hartwig.hmftools.common.sv.linx.LinxCluster.DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.Collections;
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
    private static final String SOURCE_SUBFIELD_ITEM_DELIMITER = "|";
    private static final String SOURCE_SUBFIELD_DELIMITER = ";";

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
                .add("eventIsHighDriver")
                .add("germline")
                .add("reported")
                .add("treatment")
                .add("onLabel")
                .add("level")
                .add("direction")
                .add("sources")
                .toString();
    }

    @NotNull
    private static String toLine(@NotNull ProtectEvidence evidence) {
        return new StringJoiner(FIELD_DELIMITER).add(nullToEmpty(evidence.gene()))
                .add(evidence.transcript())
                .add(String.valueOf(evidence.isCanonical()))
                .add(evidence.event())
                .add(nullToEmpty(evidence.eventIsHighDriver()))
                .add(String.valueOf(evidence.germline()))
                .add(String.valueOf(evidence.reported()))
                .add(evidence.treatment())
                .add(String.valueOf(evidence.onLabel()))
                .add(evidence.level().toString())
                .add(evidence.direction().toString())
                .add(toData(evidence.protectSources()))
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

    public static Set<ProtectSource> fromData(final String data) {
        Set<ProtectSource> protectSources = Sets.newHashSet();
        String[] dataEntries = data.split(SOURCE_SUBFIELD_DELIMITER);
        for (String entry : dataEntries) {
            String[] items = entry.split("\\" + SOURCE_SUBFIELD_ITEM_DELIMITER, -1);
            String[] urls = items[2].split(SUBFIELD_DELIMITER);
            Set<String> urlSet = Sets.newHashSet();

            for (String url : urls) {
                urlSet.add(url);
            }

            String[] urlsEvidence = items[5].split(SUBFIELD_DELIMITER);
            Set<String> urlSetEvidence = Sets.newHashSet();
            for (String url : urlsEvidence) {
                urlSetEvidence.add(url);
            }

            protectSources.add(ImmutableProtectSource.builder()
                    .source(Knowledgebase.lookupKnowledgebase(items[0]))
                    .sourceEvent(items[1])
                    .sourceUrls(urlSet.stream().sorted().collect(Collectors.toList()))
                    .evidenceType(ProtectEvidenceType.valueOf(items[3]))
                    .rangeRank(emptyToNullInteger(items[4] == null ? Strings.EMPTY : items[4]))
                    .evidenceUrls(urlSetEvidence)
                    .build());
        }
        return protectSources;
    }

    public static String toData(final Set<ProtectSource> protectSources) {
        StringBuilder sb = new StringBuilder();
        for (ProtectSource source : protectSources) {

            StringJoiner urls = new StringJoiner(SUBFIELD_DELIMITER);
            StringJoiner evidenceUrls = new StringJoiner(SUBFIELD_DELIMITER);
            sb.append(source.source().technicalDisplay()).append(SOURCE_SUBFIELD_ITEM_DELIMITER);
            sb.append(source.sourceEvent()).append(SOURCE_SUBFIELD_ITEM_DELIMITER);
            for (String sourceUrl : source.sourceUrls()) {
                urls.add(sourceUrl);
            }
            sb.append(urls).append(SOURCE_SUBFIELD_ITEM_DELIMITER);
            sb.append(source.evidenceType()).append(SOURCE_SUBFIELD_ITEM_DELIMITER);
            sb.append(source.rangeRank() == null ? Strings.EMPTY : String.valueOf(source.rangeRank()))
                    .append(SOURCE_SUBFIELD_ITEM_DELIMITER);

            for (String url : source.evidenceUrls()) {
                evidenceUrls.add(url);
            }

            sb.append(evidenceUrls).append(SOURCE_SUBFIELD_ITEM_DELIMITER);
            sb.append(SOURCE_SUBFIELD_DELIMITER);
        }

        return sb.toString();
    }

    @NotNull
    private static ProtectEvidence fromLine(@NotNull Map<String, Integer> fields, @NotNull String line) {
        String[] values = line.split(FIELD_DELIMITER, -1);

        String eventIsHighDriverField = values[fields.get("eventIsHighDriver")];
        Boolean eventIsHighDriver = !eventIsHighDriverField.isEmpty() ? Boolean.parseBoolean(eventIsHighDriverField) : null;

        String eventIsCanonicaField = values[fields.get("isCanonical")];
        Boolean isCanonical = !eventIsCanonicaField.isEmpty() ? Boolean.parseBoolean(eventIsCanonicaField) : null;

        String transcriptField = values[fields.get("transcript")];
        String transcript = !transcriptField.isEmpty() ? transcriptField : Strings.EMPTY;

        Set<ProtectSource> sources = fromData(values[fields.get("sources")]);
        return ImmutableProtectEvidence.builder()
                .gene(emptyToNullString(values[fields.get("gene")]))
                .transcript(emptyToNullString(transcript))
                .isCanonical(isCanonical)
                .event(values[fields.get("event")])
                .eventIsHighDriver(eventIsHighDriver)
                .germline(Boolean.parseBoolean(values[fields.get("germline")]))
                .reported(Boolean.parseBoolean(values[fields.get("reported")]))
                .treatment(values[fields.get("treatment")])
                .onLabel(Boolean.parseBoolean(values[fields.get("onLabel")]))
                .level(EvidenceLevel.valueOf(values[fields.get("level")]))
                .direction(EvidenceDirection.valueOf(values[fields.get("direction")]))
                .protectSources(sources)
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