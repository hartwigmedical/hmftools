package com.hartwig.hmftools.purple.qc;

import com.hartwig.hmftools.common.purple.gender.Gender;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleQC {

    private static final int SEGMENT_THRESHOLD = 150;

    public boolean overallPass() {
        return segmentPass() && genderPass();
    }

    public boolean segmentPass() {
        return segmentScore() <= SEGMENT_THRESHOLD;
    }

    public boolean genderPass() {
        return cobaltGender().equals(purpleGender());
    }

    public int segmentScore() {
        return (int) Math.round(Math.pow(ratioSegments(), 3) / Math.pow(trailingSegments(), 2) / ploidy());
    }

    public abstract int ratioSegments();

    public abstract int trailingSegments();

    public abstract double ploidy();

    public abstract Gender cobaltGender();

    public abstract Gender purpleGender();
}
