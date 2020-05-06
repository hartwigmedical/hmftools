package com.hartwig.hmftools.knowledgebasegenerator.actionability.signature;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class SignatureEvidenceAnalyzerFactory {

    private static final String DELIMITER = "\t";

    private SignatureEvidenceAnalyzerFactory(){

    }

    @NotNull
    public static SignatureEvidenceAnalyzer loadFromFileSignature(@NotNull String actionableSignatureTsv) throws IOException {
        final List<ActionableSignature> signatures = Lists.newArrayList();
        final List<String> lineSignature = Files.readAllLines(new File(actionableSignatureTsv).toPath());

        // Skip header line for signatures
        for (String lineSignatures : lineSignature.subList(1, lineSignature.size())) {
            signatures.add(fromLineSignatures(lineSignatures));
        }
        return new SignatureEvidenceAnalyzer(signatures);
    }

    @NotNull
    private static ActionableSignature fromLineSignatures(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionableSignature.builder()
                .signature(values[0])
                .source(values[1])
                .drug(values[2])
                .drugType(values[3])
                .cancerType(values[4])
                .level(values[5])
                .direction(values[6])
                .link(values[7])
                .build();
    }
}
