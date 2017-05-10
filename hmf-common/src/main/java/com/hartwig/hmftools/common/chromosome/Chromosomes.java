package com.hartwig.hmftools.common.chromosome;

import org.jetbrains.annotations.NotNull;

public final class Chromosomes {

    private Chromosomes() {
    }

    public static int asInt(@NotNull final String chromosome) {
        switch (chromosome) {
            case "X": return 23;
            case "Y": return 24;
            case "MT": return 25;
        }
        return Integer.valueOf(chromosome);
    }

    public static long length(@NotNull final String chromosome) {
        switch (chromosome) {
            case "1": return 249250622;
            case "10": return 135534749;
            case "11": return 135006518;
            case "12": return 133851897;
            case "13": return 115169880;
            case "14": return 107349541;
            case "15": return 102531394;
            case "16": return 90354755;
            case "17": return 81195212;
            case "18": return 78077250;
            case "19": return 59128985;
            case "2": return 243199374;
            case "20": return 63025522;
            case "21": return 48129897;
            case "22": return 51304568;
            case "3": return 198022431;
            case "4": return 191154277;
            case "5": return 180915261;
            case "6": return 171115068;
            case "7": return 159138664;
            case "8": return 146364023;
            case "9": return 141213432;
            case "MT": return 16571;
            case "X": return 155270561;
            case "Y": return 59373567;
        }
        return 0;
    }
}
