package com.hartwig.hmftools.datamodel.finding;

import java.util.List;

import com.hartwig.hmftools.datamodel.purple.HotspotType;
import com.hartwig.hmftools.datamodel.purple.PurpleAllelicDepth;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleGenotypeStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface SmallVariant extends Driver
{
    @NotNull
    PurpleVariant purpleVariant();

    @Nullable PurpleDriver driver();

    @NotNull
    PurpleTranscriptImpact transcriptImpact();

    @Value.Derived
    @NotNull
    default String gene()
    {
        return purpleVariant().gene();
    }

    @Value.Derived
    @NotNull
    default String chromosome()
    {
        return purpleVariant().chromosome();
    }

    @Value.Derived
    default int position()
    {
        return purpleVariant().position();
    }

    @Value.Derived
    default boolean reported()
    {
        return purpleVariant().reported();
    }

    boolean isCanonical();

    @Value.Derived
    @Nullable
    default Integer affectedCodon()
    {
        return transcriptImpact().affectedCodon();
    }

    @Value.Derived
    default double variantCopyNumber()
    {
        PurpleVariant v = purpleVariant();
        return v.adjustedCopyNumber() * Math.max(0, Math.min(1, v.adjustedVAF()));
    }

    @Value.Derived
    default double totalCopyNumber()
    {
        return purpleVariant().adjustedCopyNumber();
    }

    @Value.Derived
    default double minorAlleleCopyNumber()
    {
        return purpleVariant().minorAlleleCopyNumber();
    }

    @Value.Derived
    default boolean biallelic()
    {
        return purpleVariant().biallelic();
    }

    @Value.Derived
    @Nullable
    default Double biallelicProbability()
    {
        return purpleVariant().biallelicProbability();
    }

    @Value.Derived
    @NotNull
    default HotspotType hotspot()
    {
        return purpleVariant().hotspot();
    }

    @Value.Derived
    @Nullable
    default Double driverLikelihood()
    {
        PurpleDriver driver = driver();
        return driver != null ? driver.driverLikelihood() : null;
    }

    @Value.Derived
    default double clonalLikelihood()
    {
        return 1 - purpleVariant().subclonalLikelihood();
    }

    @Value.Derived
    @Nullable
    default List<Integer> localPhaseSets()
    {
        return purpleVariant().localPhaseSets();
    }

    @Value.Derived
    @NotNull
    default PurpleAllelicDepth tumorDepth()
    {
        return purpleVariant().tumorDepth();
    }

    @Value.Derived
    @Nullable
    default PurpleAllelicDepth rnaDepth()
    {
        return purpleVariant().rnaDepth();
    }

    @Value.Derived
    @NotNull
    default PurpleGenotypeStatus genotypeStatus()
    {
        return purpleVariant().genotypeStatus();
    }
}
