package com.hartwig.hmftools.common.variant;

public interface PurityAdjustedSomaticVariant extends Variant {

    double adjustedCopyNumber();

    double adjustedVAF();

    Clonality clonality();

    boolean lossOfHeterozygosity();
}
