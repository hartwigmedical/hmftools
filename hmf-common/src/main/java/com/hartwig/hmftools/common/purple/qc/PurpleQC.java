package com.hartwig.hmftools.common.purple.qc;

import java.util.Set;

import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.purple.gender.Gender;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleQC {

    private static final int SEGMENT_THRESHOLD = 220;
    private static final int DELETED_GENES_THRESHOLD = 280;

    @NotNull
    public PurpleQCStatus status() {
        if (!segmentPass()) {
            return PurpleQCStatus.FAIL_SEGMENT;
        }

        if (!genderPass()) {
            return PurpleQCStatus.FAIL_GENDER;
        }

        if (!deletedGenesPass()) {
            return PurpleQCStatus.FAIL_DELETED_GENES;
        }

        return PurpleQCStatus.PASS;
    }

    boolean segmentPass() {
        return segmentScore() <= SEGMENT_THRESHOLD;
    }

    boolean genderPass() {
        return cobaltGender().equals(amberGender()) || germlineAberrations().contains(GermlineAberration.KLINEFELTER);
    }

    boolean deletedGenesPass() {
        return deletedGenes() <= DELETED_GENES_THRESHOLD;
    }

    public int segmentScore() {
        return unsupportedSegments();
    }

    abstract int unsupportedSegments();

    public abstract int deletedGenes();

    public abstract double ploidy();

    @NotNull
    public abstract Gender cobaltGender();

    @NotNull
    public abstract Gender amberGender();

    @NotNull
    public abstract Set<GermlineAberration> germlineAberrations();
}
