package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface SomaticVariant extends VariantProxy
{
    boolean recovered();

    @NotNull
    String kataegis();

    double subclonalLikelihood();

    default double clonalLikelihood()
    {
        return 1 - subclonalLikelihood();
    }

    @Nullable
    AllelicDepth referenceDepth();

    double gnomadFrequency();
    SomaticLikelihood somaticLikelihood();
}