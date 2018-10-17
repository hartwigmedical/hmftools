package com.hartwig.hmftools.svanalysis.svgraph;

public enum SimplificationType {
    SimpleDeletion,
    SimpleDuplication,
    //Inversion,
    //BalancedDoubleStrandedBreakTranslocation,
    //SimpleInsertion,
    //TranslocationInsertion,
    /**
     * Adjacent events which are linked together on the same chromatid
     */
    Chain,
    /**
     * Event is linked to both sides of a fold-back inversion
     */
    ChainToFoldBackInversion,
}