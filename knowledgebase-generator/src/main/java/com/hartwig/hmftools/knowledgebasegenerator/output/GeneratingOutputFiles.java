package com.hartwig.hmftools.knowledgebasegenerator.output;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.knowledgebasegenerator.AllGenomicEvents;
import com.hartwig.hmftools.knowledgebasegenerator.cnv.KnownAmplificationDeletion;
import com.hartwig.hmftools.knowledgebasegenerator.eventtype.EventType;
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
    private static final String KNOWN_FUSION_PAIRS_INFO_TSV = "knownFusionPairsInfo.tsv";
    private static final String UNIQUE_KNOWN_FUSION_PAIRS_TSV = "uniqueKnownFusionPairs.tsv";
    private static final String UNIQUE_KNOWN_FUSION_PROMISCUOUS_THREE_TSV = "uniqueKnownFusionPromiscuousThree.tsv";
    private static final String UNIQUE_KNOWN_FUSION_PROMISCUOUS_FIVE_TSV = "uniqueKnownFusionPromiscuousFive.tsv";
    private static final String EVENT_TYPES_TSV = "eventTypes.tsv";

    public static void generatingOutputFiles(@NotNull String outputDir, @NotNull AllGenomicEvents genomicEvents) throws IOException {
        generateEventTypes(outputDir + File.separator + EVENT_TYPES_TSV, genomicEvents);
        generateUniqueKnownAmplification(outputDir + File.separator + UNIQUE_KNOWN_AMPLIFICATION_TSV, genomicEvents);
        generateUniqueKnownDeletions(outputDir + File.separator + UNIQUE_KNOWN_DELETION_TSV, genomicEvents);
        generateInfoKnownAmplification(outputDir + File.separator + KNOWN_AMPLIFICATION_INFO_TSV, genomicEvents);
        generateInfoKnownDeletions(outputDir + File.separator + KNOWN_DELETION_INFO_TSV, genomicEvents);
        generateActionableCNV(outputDir + File.separator + ACTIONABLE_CNV_TSV, genomicEvents);
        generateSignatures(outputDir + File.separator + SIGNATURES_TSV, genomicEvents);
        generateInfoKnownFusionPairs(outputDir + File.separator + KNOWN_FUSION_PAIRS_INFO_TSV, genomicEvents);
        genertateUniqueKnownFusionPairs(outputDir + File.separator + UNIQUE_KNOWN_FUSION_PAIRS_TSV, genomicEvents);
        genertateUniqueKnownFusionPromiscuousThree(outputDir + File.separator + UNIQUE_KNOWN_FUSION_PROMISCUOUS_THREE_TSV, genomicEvents);
        genertateUniqueKnownFusionPromiscuousFive(outputDir + File.separator + UNIQUE_KNOWN_FUSION_PROMISCUOUS_FIVE_TSV, genomicEvents);

    }

    private static void generateEventTypes(@NotNull String outputFile, @NotNull AllGenomicEvents genomicEvents) throws IOException {
        String headerEvents =
                "EventMap" + DELIMITER + "Map gene" + DELIMITER + "Map event" + DELIMITER
                        + "Source" + NEW_LINE;

        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, true));
        writer.write(headerEvents);

        for (EventType eventType : genomicEvents.eventType()) {
            for (Map.Entry<String, List<String>> entryDB : eventType.eventMap().entrySet()) {
                for (String event : entryDB.getValue()) {
                    writer.write(eventType.eventMap() + DELIMITER
                            + entryDB.getKey() + DELIMITER + event + DELIMITER + eventType.source() + NEW_LINE);
                }
            }
        }
        writer.close();
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

    private static void generateInfoKnownFusionPairs(@NotNull String outputFile, @NotNull AllGenomicEvents genomicEvents)
            throws IOException {
        String headerknownCNV = "gene" + DELIMITER + "eventType" + DELIMITER + "Source" + DELIMITER + "Link" + NEW_LINE;

        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, true));
        writer.write(headerknownCNV);
        for (KnownFusions knownFusionsPairs : genomicEvents.knownFusionPairs()) {
            writer.write(knownFusionsPairs.gene() + DELIMITER + knownFusionsPairs.eventType() + DELIMITER + knownFusionsPairs.source()
                    + DELIMITER + knownFusionsPairs.sourceLink() + NEW_LINE);
        }
        writer.close();
    }

    private static void genertateUniqueKnownFusionPairs(@NotNull String outputFile, @NotNull AllGenomicEvents genomicEvents)
            throws IOException {
        String headerknownFusionPairs = "Gene" + NEW_LINE;
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, true));
        writer.write(headerknownFusionPairs);
        for (String knownFusions : genomicEvents.uniqueKnownFusionPairs()) {
            writer.write(knownFusions + NEW_LINE);
        }
        writer.close();
    }

    private static void genertateUniqueKnownFusionPromiscuousThree(@NotNull String outputFile, @NotNull AllGenomicEvents genomicEvents)
            throws IOException {
        String headerknownFusionPromiscuousThree = "Gene" + NEW_LINE;
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, true));
        writer.write(headerknownFusionPromiscuousThree);
        for (String promiscuousThree : genomicEvents.knownFusionPromiscuousThree()) {
            writer.write(promiscuousThree + NEW_LINE);
        }
        writer.close();
    }

    private static void genertateUniqueKnownFusionPromiscuousFive(@NotNull String outputFile, @NotNull AllGenomicEvents genomicEvents)
            throws IOException {
        String headerknownFusionPromiscuousFive = "Gene" + NEW_LINE;
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, true));
        writer.write(headerknownFusionPromiscuousFive);
        for (String promiscuousFive : genomicEvents.knownFusionPromiscuousFive()) {
            writer.write(promiscuousFive + NEW_LINE);
        }
        writer.close();
    }
}
