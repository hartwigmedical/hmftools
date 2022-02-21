package com.hartwig.hmftools.serve.extraction.events;

import static com.hartwig.hmftools.serve.actionability.util.ActionableFileFunctions.FIELD_DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class EventInterpretationFile {

    private static final String EVENT_INTERPRETATION_TSV = "EventInterpretation.tsv";

    private EventInterpretationFile() {
    }

    @NotNull
    public static String eventInterpretationTsv(@NotNull String serveActionabilityDir) {
        return serveActionabilityDir + File.separator + EVENT_INTERPRETATION_TSV;
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
        return new StringJoiner(FIELD_DELIMITER).add("Source")
                .add("sourceEvent")
                .add("interpretGene")
                .add("interpretEvent")
                .add("interpretEventType")
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
        return new StringJoiner(FIELD_DELIMITER).add(eventInterpretation.knowledgebase().toString())
                .add(eventInterpretation.sourceEvent())
                .add(eventInterpretation.interpretGene())
                .add(eventInterpretation.interpretEvent())
                .add(eventInterpretation.interpretEventType().toString())
                .toString();
    }
}