package com.hartwig.hmftools.datamodel.immuno;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface ImmuneEscapeRecord
{
    boolean hasHlaEscape();

    boolean hasAntigenPresentationPathwayEscape();

    boolean hasIFNGammaPathwayEscape();

    boolean hasPDL1OverexpressionEscape();

    boolean hasCD58InactivationEscape();

    boolean hasEpigeneticSETDB1Escape();
}
