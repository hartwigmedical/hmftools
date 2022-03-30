package com.hartwig.hmftools.common.purple;

import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.purple.purity.FittedPurityMethod;
import com.hartwig.hmftools.common.utils.Doubles;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleQC {

    public static final int MAX_UNSUPPORTED_SEGMENTS = 220;
    public static final int MAX_DELETED_GENES = 280;
    public static final double MAX_CONTAMINATION = 0.1;
    public static final double MIN_PURITY = 0.2;

    @NotNull
    public Set<PurpleQCStatus> status() {
        Set<PurpleQCStatus> result = Sets.newHashSet();

        if (!genderPass()) {
            result.add(PurpleQCStatus.WARN_GENDER_MISMATCH);
        }

        if (unsupportedCopyNumberSegments() > MAX_UNSUPPORTED_SEGMENTS) {
            result.add(PurpleQCStatus.WARN_HIGH_COPY_NUMBER_NOISE);
        }

        if (deletedGenes() > MAX_DELETED_GENES) {
            result.add(PurpleQCStatus.WARN_DELETED_GENES);
        }

        if (Doubles.lessThan(purity(), MIN_PURITY) && !method().equals(FittedPurityMethod.NO_TUMOR)) {
            result.add(PurpleQCStatus.WARN_LOW_PURITY);
        }

        if (Doubles.greaterThan(contamination(), MAX_CONTAMINATION)) {
            result.add(PurpleQCStatus.FAIL_CONTAMINATION);
        }

        if (method().equals(FittedPurityMethod.NO_TUMOR)) {
            result.add(PurpleQCStatus.FAIL_NO_TUMOR);
        }

        if (result.isEmpty()) {
            result.add(PurpleQCStatus.PASS);
        }

        return result;
    }

    public abstract FittedPurityMethod method();

    public boolean pass() {
        return status().size() == 1 && status().contains(PurpleQCStatus.PASS);
    }

    boolean genderPass()
    {
        return cobaltGender().equals(amberGender()) || germlineAberrations().contains(GermlineAberration.KLINEFELTER);
    }

    public abstract int copyNumberSegments();

    public abstract int unsupportedCopyNumberSegments();

    public abstract int deletedGenes();

    public abstract double purity();

    public abstract double contamination();

    public abstract int amberMeanDepth();

    @NotNull
    public abstract Gender cobaltGender();

    @NotNull
    public abstract Gender amberGender();

    @NotNull
    public abstract Set<GermlineAberration> germlineAberrations();

    @Override
    public String toString() {
        return status().stream().map(Enum::toString).collect(Collectors.joining(","));
    }
}
