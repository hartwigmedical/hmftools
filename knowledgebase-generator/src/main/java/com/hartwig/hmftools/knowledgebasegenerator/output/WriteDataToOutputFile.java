package com.hartwig.hmftools.knowledgebasegenerator.output;

import java.io.BufferedWriter;
import java.io.IOException;

import com.hartwig.hmftools.knowledgebasegenerator.cnv.KnownAmplificationDeletion;

import org.jetbrains.annotations.NotNull;

public class WriteDataToOutputFile {

    private static final String DELIMITER = "\t";
    private static final String SOURCE_LINK_SEPARATOR = ";";
    private static final String NEW_LINE = "\n";

    public static void writeKnownAmplifications(@NotNull KnownAmplificationDeletion knownAmplification,
            @NotNull BufferedWriter writerKnownAmplification) throws IOException {
        writerKnownAmplification.write(
                knownAmplification.gene() + DELIMITER + knownAmplification.eventType() + DELIMITER + knownAmplification.source()
                        + NEW_LINE);
    }

    public static void writeKnownDeletion(@NotNull KnownAmplificationDeletion knownDeletion, @NotNull BufferedWriter writerKnownDeletion)
            throws IOException {
        writerKnownDeletion.write(
                knownDeletion.gene() + DELIMITER + knownDeletion.eventType() + DELIMITER + knownDeletion.source() + NEW_LINE);
    }
}
