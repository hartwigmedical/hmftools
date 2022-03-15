package com.hartwig.hmftools.serve.cancertype;

import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

public final class CancerTypeFactory {

    private static final String DELIMITER_MULTIPLE_TYPES = ";";
    private static final String SEPARATOR_NAME_AND_DOID = ",";

    private CancerTypeFactory() {
    }

    @NotNull
    public static String toString(@NotNull Set<CancerType> cancerTypes) {
        StringJoiner joiner = new StringJoiner(DELIMITER_MULTIPLE_TYPES);
        for (CancerType cancerType : cancerTypes) {
            joiner.add(cancerType.name() + SEPARATOR_NAME_AND_DOID + cancerType.doid());
        }
        return joiner.toString();
    }

    @NotNull
    public static Set<CancerType> fromString(@NotNull String cancerTypeString) {
        Set<CancerType> cancerTypes = Sets.newConcurrentHashSet();
        if (!cancerTypeString.isEmpty()) {
            String[] splitCancerTypeString = cancerTypeString.split(DELIMITER_MULTIPLE_TYPES);

            for (String cancerTypeEntry : splitCancerTypeString) {
                String[] nameAndDoid = cancerTypeEntry.split(SEPARATOR_NAME_AND_DOID);
                cancerTypes.add(ImmutableCancerType.builder().name(nameAndDoid[0]).doid(nameAndDoid[1]).build());
            }
        }

        return cancerTypes;
    }
}