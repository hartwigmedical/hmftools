package com.hartwig.hmftools.datamodel.finding;

import java.util.List;

import com.hartwig.hmftools.datamodel.driver.Driver;
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

    @NotNull PurpleDriver driver();

    @NotNull
    PurpleTranscriptImpact transcriptImpact();

    @NotNull
    PurpleTranscriptImpact otherImpact();

    @NotNull
    default String gene()
    {
        return purpleVariant().gene();
    }

    @NotNull
    default String chromosome()
    {
        return purpleVariant().chromosome();
    }

    default int position()
    {
        return purpleVariant().position();
    }

    default boolean reported()
    {
        return purpleVariant().reported();
    }

    boolean isCanonical();

    @Nullable
    default Integer affectedCodon()
    {
        return transcriptImpact().affectedCodon();
    }

    default double variantCopyNumber()
    {
        PurpleVariant v = purpleVariant();
        return v.adjustedCopyNumber() * Math.max(0, Math.min(1, v.adjustedVAF()));
    }

    default double totalCopyNumber()
    {
        return purpleVariant().adjustedCopyNumber();
    }

    default double minorAlleleCopyNumber()
    {
        return purpleVariant().minorAlleleCopyNumber();
    }

    default boolean biallelic()
    {
        return purpleVariant().biallelic();
    }

    @Nullable
    default Double biallelicProbability()
    {
        return purpleVariant().biallelicProbability();
    }

    @NotNull
    default HotspotType hotspot()
    {
        return purpleVariant().hotspot();
    }

    @NotNull
    default Double driverLikelihood()
    {
        return driver().driverLikelihood();
    }

    @Value.Derived
    default double clonalLikelihood()
    {
        return 1 - purpleVariant().subclonalLikelihood();
    }

    @Nullable
    default List<Integer> localPhaseSets()
    {
        return purpleVariant().localPhaseSets();
    }

    @NotNull
    default PurpleAllelicDepth tumorDepth()
    {
        return purpleVariant().tumorDepth();
    }

    @Nullable
    default PurpleAllelicDepth rnaDepth()
    {
        return purpleVariant().rnaDepth();
    }

    @NotNull
    default PurpleGenotypeStatus genotypeStatus()
    {
        return purpleVariant().genotypeStatus();
    }
}
