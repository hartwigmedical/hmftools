package com.hartwig.hmftools.patientreporter.cfreport.data;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class GeneUtil {

    private GeneUtil() {
    }

    @NotNull
    public static String ploidyToCopiesString(@Nullable Double ploidy, boolean hasReliablePurityFit) {
        if (!hasReliablePurityFit) {
            return DataUtil.NA_STRING;
        } else {
            return ploidy != null ? ReportResources.decimalFormat("#.#").format(ploidy) : Strings.EMPTY;
        }
    }

    @NotNull
    static String zeroPrefixed(@NotNull String location) {
        // First remove q or p arm if present.
        int armStart = location.indexOf("q");
        if (armStart < 0) {
            armStart = location.indexOf("p");
        }

        String chromosome = armStart > 0 ? location.substring(0, armStart) : location;

        try {
            int chromosomeIndex = Integer.parseInt(chromosome);
            if (chromosomeIndex < 10) {
                return "0" + location;
            } else {
                return location;
            }
        } catch (NumberFormatException exception) {
            // No need to prefix Y/X chromosomes
            return location;
        }
    }
}
