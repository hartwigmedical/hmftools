package com.hartwig.hmftools.purple.sv;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

class SampleNames {

    private final String normal;
    private final String tumor;

    SampleNames(final Set<String> names) {
        Preconditions.checkArgument(names.size() == 2);
        List<String> nameList = Lists.newArrayList(names);
        String[] nameArray = sortReferenceAndTumor(nameList.get(0), nameList.get(1));
        normal = nameArray[0];
        tumor = nameArray[1];
    }

    public String normal() {
        return normal;
    }

    public String tumor() {
        return tumor;
    }

    @VisibleForTesting
    static String[] sortReferenceAndTumor(@NotNull final String first, @NotNull final String second) {
        int intersectionLength = sampleIntersectionLength(first, second);
        if (intersectionLength != 0 && second.length() > intersectionLength) {
            final String secondRemainder = second.substring(intersectionLength);
            if (isBloodOrReference(secondRemainder)) {
                return new String[] { second, first };
            }
        }

        return new String[] { first, second };
    }

    private static boolean isBloodOrReference(String header) {
        return header.toLowerCase().startsWith("bl") || header.toLowerCase().startsWith("r");
    }

    private static int sampleIntersectionLength(@NotNull final String first, @NotNull final String second) {
        int maxIndex = Math.min(first.length(), second.length());

        for (int i = 0; i < maxIndex; i++) {
            if (first.charAt(i) != second.charAt(i)) {
                return i;
            }
        }

        return maxIndex;
    }


}
