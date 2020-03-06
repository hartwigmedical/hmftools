package com.hartwig.hmftools.knowledgebasegenerator.output;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.knowledgebasegenerator.AllGenomicEvents;
import com.hartwig.hmftools.knowledgebasegenerator.GenomicEvents;
import com.hartwig.hmftools.knowledgebasegenerator.cnv.KnownAmplificationDeletion;

import org.jetbrains.annotations.NotNull;

public class GeneratingOutputFiles {

    private static final String DELIMITER = "\t";
    private static final String NEW_LINE = "\n";

    private static final String KNOWN_AMPLIFICATION_TSV = "knownAmplification.tsv";
    private static final String KNOWN_DELETION_TSV = "knownDeletion.tsv";
    private static final String ACTIONABLE_CNV_TSV = "actionableCNV.tsv";

    public static void generatingOutputFiles(@NotNull String outputDir, @NotNull AllGenomicEvents genomicEvents) throws IOException {
        generateKnownAmplification(outputDir + File.separator + KNOWN_AMPLIFICATION_TSV, genomicEvents);
        generateKnownDeletions(outputDir + File.separator + KNOWN_DELETION_TSV, genomicEvents);
        generateActionableCNV(outputDir + File.separator + ACTIONABLE_CNV_TSV, genomicEvents);
    }

    private static void generateKnownAmplification(@NotNull String outputFile, @NotNull AllGenomicEvents genomicEvents) throws IOException {
        String headerknownCNV = "Gene" + NEW_LINE;

        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, true));
        writer.write(headerknownCNV);
        for (String amplification : genomicEvents.uniqueAmplification()) {
            writer.write(amplification + NEW_LINE);
        }
        writer.close();
    }

    private static void generateKnownDeletions(@NotNull String outputFile, @NotNull AllGenomicEvents genomicEvents) throws IOException {
        String headerknownCNV = "Gene" + NEW_LINE;
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, true));
        writer.write(headerknownCNV);
        for (String deletion : genomicEvents.uniqueDeletions()) {
            writer.write(deletion + NEW_LINE);
        }
        writer.close();
    }

    private static void generateActionableCNV(@NotNull String outputFile, @NotNull AllGenomicEvents genomicEvents) throws IOException {
        String headerActionableCNV =
                "Gene" + DELIMITER + "Type" + DELIMITER + "Source" + DELIMITER + "Links" + DELIMITER + "Drug" + DELIMITER + "Drug Type"
                        + DELIMITER + "Cancer Type" + DELIMITER + "Level" + DELIMITER + "Direction" + NEW_LINE;
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, true));
        writer.write(headerActionableCNV);
        //TODO determine actionable CNVs
        writer.write("");
        writer.close();

    }
}
