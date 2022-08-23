package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.GeneCopyNumber;

import org.jetbrains.annotations.NotNull;

public class LohGenes {

    public static long round(double copyNumber) {
        return Math.round(Math.max(0, copyNumber));
    }

    @NotNull
    public static List<GeneCopyNumber> sort(@NotNull List<GeneCopyNumber> lohGenes) {
        return lohGenes.stream().sorted((geneCopyNumber1, geneCopyNumber2) -> {
            String location1 = zeroPrefixed(geneCopyNumber1.chromosome() + geneCopyNumber1.chromosomeBand());
            String location2 = zeroPrefixed(geneCopyNumber2.chromosome() + geneCopyNumber2.chromosomeBand());

            if (location1.equals(location2)) {
                return geneCopyNumber1.geneName().compareTo(geneCopyNumber2.geneName());
            } else {
                return location1.compareTo(location2);
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    public static String zeroPrefixed(@NotNull String location) {
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
