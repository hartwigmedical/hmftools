package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;

import java.util.List;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.purple.GermlineStatus;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface Variant extends GenomePosition
{
    @NotNull
    VariantType type();

    @NotNull
    String gene();

    @NotNull
    String ref();

    @NotNull
    String alt();

    @NotNull
    AllelicDepth allelicDepth();

    @NotNull
    String canonicalTranscript();

    @NotNull
    String canonicalEffect();

    @NotNull
    CodingEffect canonicalCodingEffect();

    @NotNull
    String canonicalHgvsCodingImpact();

    @NotNull
    String canonicalHgvsProteinImpact();

    double qual();

    double mappability();

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

    default boolean isFiltered()
    {
        return !filter().equals(PASS_FILTER);
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
    List<String> reportableTranscripts();

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

    @Nullable
    String clinvarInfo();
}
