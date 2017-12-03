package com.hartwig.hmftools.common.purple.qc;

import com.hartwig.hmftools.common.purple.gender.Gender;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleQC {

    private static final int SEGMENT_THRESHOLD = 120;

    @NotNull
    public PurpleQCStatus status() {
        if (!segmentPass()) {
            return PurpleQCStatus.FAIL_SEGMENT;
        }

        if (!genderPass()) {
            return PurpleQCStatus.FAIL_GENDER;
        }

        return PurpleQCStatus.PASS;
    }

    boolean segmentPass() {
        return segmentScore() <= SEGMENT_THRESHOLD;
    }

    boolean genderPass() {
        return cobaltGender().equals(amberGender());
    }

    public int segmentScore() {
        return (int) Math.round(unsupportedSegments() / ploidy());
    }

    public abstract int unsupportedSegments();

    public abstract double ploidy();

    public abstract Gender cobaltGender();

    public abstract Gender amberGender();
}
