package com.hartwig.hmftools.knowledgebasegenerator.output;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.jetbrains.annotations.NotNull;

public class GeneratingOutputFiles {

    private static final String DELIMITER = "\t";
    private static final String SOURCE_LINK_SEPARATOR = ";";
    private static final String NEW_LINE = "\n";

    private static final String KNOWN_AMPLIFICATION_TSV = "knownAmplification.tsv";
    private static final String KNOWN_DELETION_TSV = "knownDeletion.tsv";
    private static final String ACTIONABLE_CNV_TSV = "actionableCNV.tsv";

    @NotNull
    public static BufferedWriter generateKnownAmplification(@NotNull String outputDir) throws IOException {
        String outputFile = outputDir + File.separator + KNOWN_AMPLIFICATION_TSV;
        String headerknownCNV = "Gene" + DELIMITER + "Type" + DELIMITER + "Source" + DELIMITER + "Links" + NEW_LINE;

        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, true));
        writer.write(headerknownCNV);
        return writer;
    }

    @NotNull
    public static BufferedWriter generateKnownDeletions(@NotNull String outputDir) throws IOException {
        String outputFile = outputDir + File.separator + KNOWN_DELETION_TSV;
        String headerknownCNV = "Gene" + DELIMITER + "Type" + DELIMITER + "Source" + DELIMITER + "Links" + NEW_LINE;
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, true));
        writer.write(headerknownCNV);
        return writer;
    }

    @NotNull
    public static BufferedWriter generateActionableCNV(@NotNull String outputDir) throws IOException{
        String outputFile = outputDir + File.separator + ACTIONABLE_CNV_TSV;
        String headerActionableCNV =
                "Gene" + DELIMITER + "Type" + DELIMITER + "Source" + DELIMITER + "Links" + DELIMITER + "Drug" + DELIMITER + "Drug Type"
                        + DELIMITER + "Cancer Type" + DELIMITER + "Level" + DELIMITER + "Direction" + NEW_LINE;
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, true));
        writer.write(headerActionableCNV);
        return writer;
    }
}
