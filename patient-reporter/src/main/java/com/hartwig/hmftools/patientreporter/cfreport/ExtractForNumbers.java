package com.hartwig.hmftools.patientreporter.cfreport;

import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ExtractForNumbers {

    private ExtractForNumbers() {
    }

    @NotNull
    public static String extractAllForNumbers(boolean isQCFail, double purity, boolean hasReliablePurity, @Nullable QCFailReason reason) {
        return isQCFail ? reason.forNumber() : successForNumber(purity, hasReliablePurity);
    }

    @NotNull
    public static String successForNumber(double purity, boolean hasReliablePurity) {
        return (hasReliablePurity && purity > 0.195) ? ForNumber.FOR_080.display() : ForNumber.FOR_103.display();
    }
}
