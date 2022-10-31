package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.patientreporter.algo.LohGenesReporting;

import org.jetbrains.annotations.NotNull;

public class LohGenes {

    @NotNull
    public static List<LohGenesReporting> sort(@NotNull List<LohGenesReporting> lohGenes) {
        return lohGenes.stream()
                .sorted(Comparator.comparing((LohGenesReporting lohGenesReporting) -> lohGenesReporting.gene()))
                .collect(Collectors.toList());
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
