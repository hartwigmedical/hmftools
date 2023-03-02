package com.hartwig.hmftools.datamodel.purple;

import com.hartwig.hmftools.datamodel.genome.chromosome.GermlineAberration;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.Set;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleQC
{

    public abstract Set<PurpleQCStatus> status();

    public abstract FittedPurityMethod method();

    public boolean pass()
    {
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
}
