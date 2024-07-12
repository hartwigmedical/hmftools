package com.hartwig.hmftools.datamodel.purple;

import java.util.Set;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface PurpleQC
{
    @NotNull
    Set<PurpleQCStatus> status();

    @NotNull
    Set<PurpleGermlineAberration> germlineAberrations();

    int amberMeanDepth();

    double contamination();

    int totalCopyNumberSegments();

    int unsupportedCopyNumberSegments();

    int deletedGenes();
}
