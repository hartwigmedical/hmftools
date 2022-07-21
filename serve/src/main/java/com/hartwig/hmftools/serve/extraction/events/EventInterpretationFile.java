package com.hartwig.hmftools.serve.extraction.events;

import static com.hartwig.hmftools.serve.actionability.util.ActionableFileFunctions.FIELD_DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.classification.EventType;

import org.jetbrains.annotations.NotNull;

public final class EventInterpretationFile {

    private static final String EVENT_INTERPRETATION_TSV = "EventInterpretation.tsv";

    private EventInterpretationFile() {
    }

    @NotNull
    public static String eventInterpretationTsv(@NotNull String serveActionabilityDir) {
        return serveActionabilityDir + File.separator + EVENT_INTERPRETATION_TSV;
    }

    @NotNull
    public static List<EventInterpretation> read(@NotNull String file) throws IOException {
        List<String> lines = Files.readAllLines(new File(file).toPath());

        return fromLines(lines.subList(1, lines.size()));
    }

    @NotNull
    @VisibleForTesting
    static List<EventInterpretation> fromLines(@NotNull List<String> lines) {
        List<EventInterpretation> eventInterpretations = Lists.newArrayList();
        for (String line : lines) {
            eventInterpretations.add(fromLine(line));
        }
        return eventInterpretations;
    }

    @NotNull
    private static EventInterpretation fromLine(@NotNull String line) {
        String[] values = line.split(FIELD_DELIMITER);

        return ImmutableEventInterpretation.builder()
                .source(Knowledgebase.lookupKnowledgebase(values[0]))
                .sourceEvent(values[1])
                .interpretedGene(values[2])
                .interpretedEvent(values[3])
                .interpretedEventType(EventType.valueOf(values[4]))
                .build();
    }

    public static void write(@NotNull String eventInterpretationTsv, @NotNull Iterable<EventInterpretation> eventInterpretations)
            throws IOException {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        lines.addAll(toLines(eventInterpretations));

        Files.write(new File(eventInterpretationTsv).toPath(), lines);
    }

    @NotNull
    private static String header() {
        return new StringJoiner(FIELD_DELIMITER).add("source")
                .add("sourceEvent")
                .add("interpretedGene")
                .add("interpretedEvent")
                .add("interpretedEventType")
                .toString();
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull Iterable<EventInterpretation> eventInterpretations) {
        List<String> lines = Lists.newArrayList();
        for (EventInterpretation eventInterpretation : eventInterpretations) {
            lines.add(toLine(eventInterpretation));
        }
        return lines;
    }

    @NotNull
    private static String toLine(@NotNull EventInterpretation eventInterpretation) {
        return new StringJoiner(FIELD_DELIMITER).add(eventInterpretation.source().toString())
                .add(eventInterpretation.sourceEvent())
                .add(eventInterpretation.interpretedGene())
                .add(eventInterpretation.interpretedEvent())
                .add(eventInterpretation.interpretedEventType().toString())
                .toString();
    }
}