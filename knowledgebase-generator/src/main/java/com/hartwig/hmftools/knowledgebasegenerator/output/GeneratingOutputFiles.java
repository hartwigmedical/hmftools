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


    public static void generatingOutputFiles(@NotNull String outputDir)  throws IOException{
        generateKnownAmplification(outputDir + File.separator + KNOWN_AMPLIFICATION_TSV);
        generateKnownDeletions(outputDir + File.separator + KNOWN_DELETION_TSV);
        generateActionableCNV(outputDir + File.separator + ACTIONABLE_CNV_TSV);
    }

    private static void generateKnownAmplification(@NotNull String outputFile) throws IOException {
        String headerknownCNV = "Gene" + DELIMITER + "Type" + DELIMITER + "Sources" + DELIMITER + "Links" + NEW_LINE;

        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, true));
        writer.write(headerknownCNV);
        writer.close();
    }

    private static void generateKnownDeletions(@NotNull String outputFile) throws IOException {
        String headerknownCNV = "Gene" + DELIMITER + "Type" + DELIMITER + "Sources" + DELIMITER + "Links" + NEW_LINE;
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, true));
        writer.write(headerknownCNV);
        writer.close();
    }

    private static void generateActionableCNV(@NotNull String outputFile) throws IOException{
        String headerActionableCNV =
                "Gene" + DELIMITER + "Type" + DELIMITER + "Sources" + DELIMITER + "Links" + DELIMITER + "Drug" + DELIMITER + "Drug Type"
                        + DELIMITER + "Cancer Type" + DELIMITER + "Level" + DELIMITER + "Direction" + NEW_LINE;
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, true));
        writer.write(headerActionableCNV);
        writer.close();

    }
}
