package com.hartwig.hmftools.orange.report.datamodel;

import java.util.List;

import com.hartwig.hmftools.datamodel.purple.HotspotType;
import com.hartwig.hmftools.datamodel.purple.PurpleAllelicDepth;
import com.hartwig.hmftools.datamodel.purple.PurpleGenotypeStatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class VariantEntry
{
    @NotNull
    public abstract String gene();

    public abstract boolean isCanonical();

    @Nullable
    public abstract Integer affectedCodon();

    @NotNull
    public abstract String impact();

    public abstract double variantCopyNumber();

    public abstract double totalCopyNumber();

    public abstract double minorAlleleCopyNumber();

    public abstract boolean biallelic();

    @Nullable
    public abstract Double biallelicProbability();

    @NotNull
    public abstract HotspotType hotspot();

    @Nullable
    public abstract Double driverLikelihood();

    public abstract double clonalLikelihood();

    @Nullable
    public abstract List<Integer> localPhaseSets();

    @Nullable
    public abstract PurpleAllelicDepth rnaDepth();

    @NotNull
    public abstract PurpleGenotypeStatus genotypeStatus();

}
