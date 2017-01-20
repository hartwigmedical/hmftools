package com.hartwig.hmftools.patientreporter.slicing;

import org.jetbrains.annotations.NotNull;

public final class HMFSlicingAnnotationTestFactory {

    private HMFSlicingAnnotationTestFactory() {
    }

    @NotNull
    public static HMFSlicingAnnotation create(@NotNull final String transcriptID, final int transcriptVersion,
            @NotNull final String gene) {
        return new HMFSlicingAnnotation(transcriptID, transcriptVersion, gene);
    }
}
