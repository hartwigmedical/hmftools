package com.hartwig.hmftools.svanalysis.svgraph;

public enum SimplificationType {
    SimpleIndel,
    SimpleDuplication,
    TranslocationInsertion,
    SimpleInversion,
    //BalancedDoubleStrandedBreakTranslocation,
    //SimpleInsertion,
    /**
     * Adjacent events which are linked together on the same chromatid
     */
    Chain,
    /**
     * Event is linked to both sides of a fold-back inversion
     */
    ChainToFoldBackInversion,
}