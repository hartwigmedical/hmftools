package com.hartwig.hmftools.common.gene;

import com.hartwig.hmftools.common.copynumber.CopyNumber;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface GeneCopyNumber extends GeneRegion, CopyNumber {

    @NotNull
    String gene();

    double maxCopyNumber();

    double minCopyNumber();

    double meanCopyNumber();

    int somaticRegions();

    int germlineHet2HomRegions();

    int germlineHomRegions();

    default int totalRegions() {
        return somaticRegions() + germlineHet2HomRegions() + germlineHomRegions();
    }

    @Override
    default int value() {
        return (int) Math.max(0, Math.round(minCopyNumber()));
    }
}
