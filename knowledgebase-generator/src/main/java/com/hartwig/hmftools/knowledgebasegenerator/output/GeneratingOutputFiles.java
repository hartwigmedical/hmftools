package com.hartwig.hmftools.knowledgebasegenerator.output;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import com.hartwig.hmftools.knowledgebasegenerator.AllGenomicEvents;
import com.hartwig.hmftools.knowledgebasegenerator.cnv.KnownAmplificationDeletion;
import com.hartwig.hmftools.knowledgebasegenerator.fusion.KnownFusions;
import com.hartwig.hmftools.knowledgebasegenerator.signatures.Signatures;

import org.jetbrains.annotations.NotNull;

public class GeneratingOutputFiles {

    private static final String DELIMITER = "\t";
    private static final String NEW_LINE = "\n";

    private static final String UNIQUE_KNOWN_AMPLIFICATION_TSV = "uniqueKnownAmplification.tsv";
    private static final String UNIQUE_KNOWN_DELETION_TSV = "uniqueKnownDeletion.tsv";
    private static final String KNOWN_AMPLIFICATION_INFO_TSV = "knownAmplificationInfo.tsv";
    private static final String KNOWN_DELETION_INFO_TSV = "knownDeletionInfo.tsv";
    private static final String ACTIONABLE_CNV_TSV = "actionableCNV.tsv";
    private static final String SIGNATURES_TSV = "signatures.tsv";
    private static final String KNOWN_FUSION_PAIRS_TSV = "knownFusionPairs.tsv";

    public static void generatingOutputFiles(@NotNull String outputDir, @NotNull AllGenomicEvents genomicEvents) throws IOException {
        generateUniqueKnownAmplification(outputDir + File.separator + UNIQUE_KNOWN_AMPLIFICATION_TSV, genomicEvents);
        generateUniqueKnownDeletions(outputDir + File.separator + UNIQUE_KNOWN_DELETION_TSV, genomicEvents);
        generateInfoKnownAmplification(outputDir + File.separator + KNOWN_AMPLIFICATION_INFO_TSV, genomicEvents);
        generateInfoKnownDeletions(outputDir + File.separator + KNOWN_DELETION_INFO_TSV, genomicEvents);
        generateActionableCNV(outputDir + File.separator + ACTIONABLE_CNV_TSV, genomicEvents);
        generateSignatures(outputDir + File.separator + SIGNATURES_TSV, genomicEvents);
        generateKnownFusionPairs(outputDir + File.separator + KNOWN_FUSION_PAIRS_TSV, genomicEvents);

    }

    private static void generateUniqueKnownAmplification(@NotNull String outputFile, @NotNull AllGenomicEvents genomicEvents)
            throws IOException {
        String headerknownCNV = "Gene" + NEW_LINE;

        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, true));
        writer.write(headerknownCNV);
        for (String amplification : genomicEvents.uniqueAmplification()) {
            writer.write(amplification + NEW_LINE);
        }
        writer.close();
    }

    private static void generateUniqueKnownDeletions(@NotNull String outputFile, @NotNull AllGenomicEvents genomicEvents)
            throws IOException {
        String headerknownCNV = "Gene" + NEW_LINE;
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, true));
        writer.write(headerknownCNV);
        for (String deletion : genomicEvents.uniqueDeletions()) {
            writer.write(deletion + NEW_LINE);
        }
        writer.close();
    }

    private static void generateInfoKnownAmplification(@NotNull String outputFile, @NotNull AllGenomicEvents genomicEvents)
            throws IOException {
        String headerknownCNV = "Gene" + DELIMITER + "Source" + NEW_LINE;

        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, true));
        writer.write(headerknownCNV);
        for (KnownAmplificationDeletion amplification : genomicEvents.knownAmplifications()) {
            writer.write(amplification.gene() + DELIMITER + amplification.source() + NEW_LINE);
        }
        writer.close();
    }

    private static void generateInfoKnownDeletions(@NotNull String outputFile, @NotNull AllGenomicEvents genomicEvents) throws IOException {
        String headerknownCNV = "Gene" + DELIMITER + "Source" + NEW_LINE;
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, true));
        writer.write(headerknownCNV);
        for (KnownAmplificationDeletion deletion : genomicEvents.knownDeletions()) {
            writer.write(deletion.gene() + DELIMITER + deletion.source() + NEW_LINE);
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

    private static void generateSignatures(@NotNull String outputFile, @NotNull AllGenomicEvents genomicEvents) throws IOException {
        String headerknownCNV = "event" + DELIMITER + "Source" + DELIMITER + "Link" + NEW_LINE;

        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, true));
        writer.write(headerknownCNV);
        for (Signatures signatures : genomicEvents.signatures()) {
            writer.write(signatures.eventType() + DELIMITER + signatures.source() + DELIMITER + signatures.sourceLink() + NEW_LINE);
        }
        writer.close();
    }

    private static void generateKnownFusionPairs(@NotNull String outputFile, @NotNull AllGenomicEvents genomicEvents) throws IOException {
        String headerknownCNV = "gene" + DELIMITER + "eventType" + DELIMITER + "Source" + DELIMITER + "Link" + NEW_LINE;

        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, true));
        writer.write(headerknownCNV);
        for (KnownFusions knownFusionsPairs : genomicEvents.knownFusionPairs()) {
            writer.write(knownFusionsPairs.gene() + DELIMITER + knownFusionsPairs.eventType() + DELIMITER + knownFusionsPairs.source()
                    + DELIMITER + knownFusionsPairs.sourceLink() + NEW_LINE);
        }
        writer.close();
    }
}
