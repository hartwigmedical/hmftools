package com.hartwig.hmftools.common.actionability.cancertype;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CancerTypeAnalyzer {
    private static final String DELIMITER = "\t";

    @NotNull
    private final List<CancerTypeReading> cancerTypeDoids;

    public CancerTypeAnalyzer(@NotNull List<CancerTypeReading> cancerTypeDoids) {
        this.cancerTypeDoids = cancerTypeDoids;
    }

    @NotNull
    public static CancerTypeAnalyzer loadFromFile(@NotNull String fileCancerType) throws IOException {
        final List<CancerTypeReading> cancerTypeWithDOID = Lists.newArrayList();
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
            return Strings.EMPTY;
        }
    }

    @VisibleForTesting
    public boolean foundTumorLocation(@NotNull String tumorLocationKnowledgebase, @Nullable String doidsPrimaryTumorLocation) {
        boolean booleanValueRange = false;
        for (CancerTypeReading cancerTypeDoidKnowledgeBase : cancerTypeDoids) {
            if (tumorLocationKnowledgebase.equals(cancerTypeDoidKnowledgeBase.cancerType())) {
                if (doidsPrimaryTumorLocation != null) {
                    if (doidsPrimaryTumorLocation.contains(";")) {
                        String[] multipleDoidsPrimaryTumorLocation = doidsPrimaryTumorLocation.split(";");
                        if (cancerTypeDoidKnowledgeBase.doidSet().contains(multipleDoidsPrimaryTumorLocation[0])
                                || cancerTypeDoidKnowledgeBase.doidSet().contains(multipleDoidsPrimaryTumorLocation[1])) {
                            booleanValueRange = true;
                        } else if (cancerTypeDoidKnowledgeBase.doidSet().contains(multipleDoidsPrimaryTumorLocation[0])
                                && cancerTypeDoidKnowledgeBase.doidSet().contains(multipleDoidsPrimaryTumorLocation[1])) {
                            booleanValueRange = true;
                        }
                    } else if (cancerTypeDoidKnowledgeBase.doidSet().contains(doidsPrimaryTumorLocation)) {
                        booleanValueRange = true;
                    }
                }
            }
        }
        return booleanValueRange;
    }
}
