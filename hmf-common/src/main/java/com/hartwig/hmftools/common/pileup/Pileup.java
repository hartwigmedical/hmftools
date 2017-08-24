package com.hartwig.hmftools.common.pileup;

import com.hartwig.hmftools.common.position.GenomePosition;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface Pileup extends GenomePosition {

    int readCount();

    @NotNull
    String referenceBase();

    int referenceCount();

    int gMismatchCount();

    int aMismatchCount();

    int tMismatchCount();

    int cMismatchCount();

    default int mismatchCount(char base) {
        switch (Character.toUpperCase(base)) {
            case 'G':
                return gMismatchCount();
            case 'A':
                return aMismatchCount();
            case 'T':
                return tMismatchCount();
            case 'C':
                return cMismatchCount();
            default:
                return 0;
        }
    }
}
