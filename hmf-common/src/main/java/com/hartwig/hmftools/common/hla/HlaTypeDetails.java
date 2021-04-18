package com.hartwig.hmftools.common.hla;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface HlaTypeDetails {

    @NotNull
    String type();

    int referenceUniqueCoverage();

    int referenceSharedCoverage();

    int referenceWildcardCoverage();

    int tumorUniqueCoverage();

    int tumorSharedCoverage();

    int tumorWildcardCoverage();

    double tumorCopyNumber();

    double somaticMissense();

    double somaticNonsenseOrFrameshift();

    double somaticSplice();

    double somaticSynonymous();

    double somaticInframeIndel();
}
