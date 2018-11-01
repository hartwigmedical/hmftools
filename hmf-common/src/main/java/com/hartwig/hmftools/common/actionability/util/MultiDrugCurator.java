package com.hartwig.hmftools.common.actionability.util;

import java.util.Arrays;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class MultiDrugCurator {

    private MultiDrugCurator() {
    }

    @NotNull
    public static String reformat(@NotNull String drug) {
        List<String> drugParts = sortedAndTrimmedParts(drug);

        StringBuilder finalDrug = new StringBuilder(drugParts.get(0));
        for (int i = 1; i < drugParts.size(); i++) {
            finalDrug.append(" + ").append(drugParts.get(i));
        }

        return finalDrug.toString();
    }

    @NotNull
    private static List<String> sortedAndTrimmedParts(@NotNull String drug) {
        List<String> drugParts = Arrays.asList(drug.split(Pattern.quote("+")));
        List<String> trimmedParts = Lists.newArrayList();
        for (String part : drugParts) {
            trimmedParts.add(part.trim());
        }
        return trimmedParts.stream().sorted().collect(Collectors.toList());
    }
}
