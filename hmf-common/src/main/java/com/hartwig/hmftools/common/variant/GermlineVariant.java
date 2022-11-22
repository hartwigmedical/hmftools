package com.hartwig.hmftools.common.variant;

import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.purple.GermlineStatus;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface GermlineVariant extends Variant
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
    VariantTier tier();

    @Nullable
    AllelicDepth rnaDepth();

    @Nullable
    String clinvarInfo();

    @Nullable
    String pathogenicity();

    boolean pathogenic();
}