package com.hartwig.hmftools.purple.somatic;

import org.apache.logging.log4j.util.Strings;

public class Microhomology {
    // length of normal should be deleted size * 2

    public static String microhomology(final String ref, final String alt, final String normal) {
        int refIndex = normal.indexOf(ref);
        String tumor = alt + normal.substring(refIndex + ref.length());
        String commonPrefix = commonPrefix(ref, tumor);
        return commonPrefix.charAt(0) == normal.charAt(ref.length() - 1) ? commonPrefix : commonPrefix.substring(1);
    }

    private static String commonPrefix(String normal, String tumor) {
        int minLength = Math.min(tumor.length(), normal.length());
        if (minLength == 0) {
            return Strings.EMPTY;
        }

        for (int i = 0; i < minLength; i++) {
            if (normal.charAt(i) != tumor.charAt(i)) {
                return i == 0 ? Strings.EMPTY : normal.substring(0, i);
            }
        }
        return normal.substring(0, minLength);
    }
}
