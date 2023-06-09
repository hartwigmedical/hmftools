package com.hartwig.hmftools.common.variant;

import java.util.List;

import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.purple.GermlineStatus;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface SomaticVariant extends Variant
{
    double qual();

    @NotNull
    String filter();

    int genesAffected();

    boolean spliceRegion();

    @NotNull
    String otherReportedEffects();

    @NotNull
    CodingEffect worstCodingEffect();

    @NotNull
    Hotspot hotspot();

    boolean recovered();

    default boolean isHotspot()
    {
        return hotspot() == Hotspot.HOTSPOT;
    }

    double mappability();

    default boolean isFiltered()
    {
        return !filter().equals(SomaticVariantFactory.PASS_FILTER);
    }

    default boolean isSnp()
    {
        return type() == VariantType.SNP;
    }

    double adjustedCopyNumber();

    double adjustedVAF();

    double minorAlleleCopyNumber();

    double variantCopyNumber();

    boolean biallelic();

    boolean reported();

    @NotNull
    GenotypeStatus genotypeStatus();

    @NotNull
    GermlineStatus germlineStatus();

    @NotNull
    String trinucleotideContext();

    @NotNull
    String microhomology();

    @NotNull
    String repeatSequence();

    int repeatCount();

    @NotNull
    String kataegis();

    @NotNull
    VariantTier tier();

    double subclonalLikelihood();

    default double clonalLikelihood()
    {
        return 1 - subclonalLikelihood();
    }

    @Nullable
    AllelicDepth rnaDepth();

    @Nullable
    AllelicDepth referenceDepth();

    @Nullable
    List<Integer> localPhaseSets();

    default Integer topLocalPhaseSet()
    {
        return localPhaseSets() != null && !localPhaseSets().isEmpty() ? localPhaseSets().get(0) : null;
    }

    default String localPhaseSetsStr()
    {
        return SomaticVariantFactory.localPhaseSetsStr(localPhaseSets());
    }

    default boolean hasLocalPhaseSets()
    {
        return localPhaseSets() != null && !localPhaseSets().isEmpty();
    }

    String clinvarInfo();
    double gnomadFrequency();
    SomaticLikelihood somaticLikelihood();
}