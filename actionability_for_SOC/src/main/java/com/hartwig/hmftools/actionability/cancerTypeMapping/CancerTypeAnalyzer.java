package com.hartwig.hmftools.actionability.cancerTypeMapping;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CancerTypeAnalyzer {
    private static final org.apache.logging.log4j.Logger LOGGER = org.apache.logging.log4j.LogManager.getLogger(CancerTypeAnalyzer.class);
    private static final String DELIMITER = "\t";

    @NotNull
    private final List<CancerTypeReading> cancerTypeDoids;

    public CancerTypeAnalyzer(@NotNull final List<CancerTypeReading> cancerTypeDoids) {
        this.cancerTypeDoids = cancerTypeDoids;
    }

    @NotNull
    public static CancerTypeAnalyzer loadFromFile(String fileCancerType) throws IOException {
        final List<CancerTypeReading> cancerTypeWithDOID = new ArrayList<>();
        final List<String> lineCancerType = Files.readAllLines(new File(fileCancerType).toPath());

        for (int i = 1; i < lineCancerType.size(); i++) {
            cancerTypeWithDOID.add(fromLine(lineCancerType.get(i)));
        }
        return new CancerTypeAnalyzer(cancerTypeWithDOID);
    }

    @NotNull
    private static CancerTypeReading fromLine(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        String doidSetValue = EmptyDoidSet(values);
        return  ImmutableCancerTypeReading.builder()
                .cancerType(values[0])
                .doidSet(doidSetValue).build();
    }

    @NotNull
    private static String EmptyDoidSet(@NotNull String[] value) {
        try {
            return value[1];
        } catch (ArrayIndexOutOfBoundsException e) {
            LOGGER.warn("IndexOutOfBoundsException: " + e.getMessage());
            return " ";
        }
    }

    @NotNull
    public boolean foundTumorLocation (@NotNull String tumorLocationKnowledgebase, @Nullable String doid) {
        Boolean booleanValueRange = false;
        for (int i = 1; i < cancerTypeDoids.size(); i ++) {
            if(cancerTypeDoids.get(i).cancerType().contains(tumorLocationKnowledgebase)){
                if(cancerTypeDoids.get(i).doidSet().contains(doid)){
                    booleanValueRange = true;
                }
            }
        }
        return booleanValueRange;
    }
}
