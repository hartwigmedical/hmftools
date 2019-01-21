package com.hartwig.hmftools.svanalysis.visualisation;

enum ChainColor {
    ;

    static String color(int chainId) {
        switch (chainId) {
            case 1:
                return "color=black";
            case 2:
                return "color=blue";
            default:
                return "color=red";
        }

    }

}
