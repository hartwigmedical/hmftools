package com.hartwig.hmftools.svanalysis.visualisation;

enum ChainColor {
    ;

    static String color(int chainId) {
        switch (chainId) {
            case 0:
                return "color=purple";
            case 1:
                return "color=red";
            case 2:
                return "color=green";
            default:
                return "color=black";
        }

    }

}
