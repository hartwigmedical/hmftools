package com.hartwig.hmftools.common.variant;

public interface PurityAdjustedSomaticVariant extends SomaticVariant {

    double adjustedCopyNumber();

    double adjustedVAF();

    Clonality clonality();

    boolean lossOfHeterozygosity();
}
