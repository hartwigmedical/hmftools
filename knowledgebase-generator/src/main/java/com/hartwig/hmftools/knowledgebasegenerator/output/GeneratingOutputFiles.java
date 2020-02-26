package com.hartwig.hmftools.knowledgebasegenerator.output;

public class GeneratingOutputFiles {

    private static final String DELIMITER = "\t";
    private static final String SOURCE_LINK_SEPARATOR = ";";
    private static final String NEW_LINE = "\n";

    public static void generatingOutputFiles() {
        generateKnownAmplification();
        generateKnownDeletions();
        generateActionableCNV();
    }

    public static void generateKnownAmplification() {
        String headerknownCNV = "Gene" + DELIMITER + "Type" + DELIMITER + "Sources" + DELIMITER + "Links" + NEW_LINE;
    }

    public static void generateKnownDeletions() {
        String headerknownCNV = "Gene" + DELIMITER + "Type" + DELIMITER + "Sources" + DELIMITER + "Links" + NEW_LINE;
    }

    public static void generateActionableCNV() {
        String headerActionableCNV =
                "Gene" + DELIMITER + "Type" + DELIMITER + "Sources" + DELIMITER + "Links" + DELIMITER + "Drug" + DELIMITER + "Drug Type"
                        + DELIMITER + "Cancer Type" + DELIMITER + "Level" + DELIMITER + "Direction" + NEW_LINE;

    }
}
