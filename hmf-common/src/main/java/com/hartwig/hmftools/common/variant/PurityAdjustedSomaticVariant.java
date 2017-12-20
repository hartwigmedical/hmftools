package com.hartwig.hmftools.common.variant;

public interface PurityAdjustedSomaticVariant extends SomaticVariant {

    double adjustedCopyNumber();

    double adjustedVAF();

    double ploidy();

    boolean lossOfHeterozygosity();
}
