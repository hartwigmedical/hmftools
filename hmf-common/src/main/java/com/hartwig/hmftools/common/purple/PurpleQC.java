package com.hartwig.hmftools.common.purple;

import static com.hartwig.hmftools.common.purple.PurpleQCStatus.STATUS_DELIM;

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

    public abstract Set<PurpleQCStatus> status();

    public abstract FittedPurityMethod method();

    public boolean pass() {
        return status().size() == 1 && status().contains(PurpleQCStatus.PASS);
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
        return status().stream().map(Enum::toString).collect(Collectors.joining(STATUS_DELIM));
    }
}
