package com.hartwig.hmftools.common.cobalt;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface CobaltRatio extends GenomePosition
{
    double referenceReadDepth();

    double tumorReadDepth();

    double referenceGCRatio();

    double referenceGCDiploidRatio();

    double tumorGCRatio();

    double referenceGcContent();

    double tumorGcContent();
}
