package com.hartwig.hmftools.actionability.cancerTypeMapping;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CancerTypeAnalyzer {
    private static final Logger LOGGER = LogManager.getLogger(CancerTypeAnalyzer.class);
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
        String doidSetValue = emptyDoidSet(values);
        return ImmutableCancerTypeReading.builder().cancerType(values[0]).doidSet(doidSetValue).build();
    }

    @NotNull
    private static String emptyDoidSet(@NotNull String[] value) {
        try {
            return value[1];
        } catch (ArrayIndexOutOfBoundsException e) {
            LOGGER.warn("IndexOutOfBoundsException: " + e.getMessage());
            return Strings.EMPTY;
        }
    }

    public boolean foundTumorLocation(@NotNull String tumorLocationKnowledgebase, @Nullable String doidsPrimaryTumorLocation) {
        Boolean booleanValueRange = false;
        for (CancerTypeReading cancerTypeDoidKnowledgeBase : cancerTypeDoids) {
            if (tumorLocationKnowledgebase.equals(cancerTypeDoidKnowledgeBase.cancerType())) {
                if (cancerTypeDoidKnowledgeBase.doidSet().contains(doidsPrimaryTumorLocation)) {
                    booleanValueRange = true;
                }
            }
        }
        return booleanValueRange;
    }
}
