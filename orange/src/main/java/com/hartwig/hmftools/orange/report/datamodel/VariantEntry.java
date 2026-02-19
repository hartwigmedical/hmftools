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
    public abstract String gene();

    public abstract int depth();
    public abstract double vaf();

    public abstract boolean isCanonical();

    @Nullable
    public abstract Integer affectedCodon();

    public abstract String impact();

    public abstract double variantCopyNumber();

    public abstract double totalCopyNumber();

    public abstract double minorAlleleCopyNumber();

    public abstract boolean biallelic();

    public abstract Double biallelicProbability();

    public abstract HotspotType hotspot();

    @Nullable
    public abstract Double driverLikelihood();

    public abstract double clonalLikelihood();

    public abstract String somaticLikelihood();

    @Nullable
    public abstract PurpleAllelicDepth rnaDepth();

    public abstract PurpleGenotypeStatus genotypeStatus();

}
