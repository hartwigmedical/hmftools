package com.hartwig.hmftools.svanalysis.visualisation;

enum ChainColor {
    ;

    static String color(int chainId) {
        switch (chainId) {
            case 1:
                return "color=blue";
            case 2:
                return "color=black";
            default:
                return "color=red";
        }

    }

}
