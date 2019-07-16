package com.hartwig.hmftools.common.variant;

import java.util.List;

import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.common.variant.cosmic.CosmicAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface SomaticVariant extends GenomePosition, AllelicDepth {

    @NotNull
    String ref();

    @NotNull
    String alt();

    @NotNull
    VariantType type();

    @NotNull
    String filter();

    @Nullable
    String dbsnpID();

    @NotNull
    List<String> cosmicIDs();

    @Nullable
    String canonicalCosmicID();

    @NotNull
    List<SnpEffAnnotation> snpEffAnnotations();

    @NotNull
    List<CosmicAnnotation> cosmicAnnotations();

    @NotNull
    String gene();

    int genesEffected();

    @NotNull
    String worstEffect();

    @NotNull
    String worstEffectTranscript();

    @NotNull
    CodingEffect worstCodingEffect();

    @NotNull
    String canonicalEffect();

    @NotNull
    CodingEffect canonicalCodingEffect();

    @NotNull
    String canonicalHgvsCodingImpact();

    @NotNull
    String canonicalHgvsProteinImpact();

    @NotNull
    Hotspot hotspot();

    boolean recovered();

    default boolean isHotspot() {
        return hotspot() == Hotspot.HOTSPOT;
    }

    double mappability();

    default boolean isDBSNP() {
        return dbsnpID() != null;
    }

    default boolean isCOSMIC() {
        return !cosmicIDs().isEmpty();
    }

    default boolean isFiltered() {
        return !filter().equals(SomaticVariantFactory.PASS_FILTER);
    }

    default boolean isSnp() {
        return type() == VariantType.SNP;
    }

    double adjustedCopyNumber();

    double adjustedVAF();

    double minorAllelePloidy();

    double ploidy();

    boolean biallelic();

    @NotNull
    GermlineStatus germlineStatus();

    @NotNull
    String trinucleotideContext();

    @NotNull
    String microhomology();

    @NotNull
    String repeatSequence();

    int repeatCount();

    boolean highConfidenceRegion();

    @NotNull
    String kataegis();

    double subclonalLikelihood();

    default double clonalLikelihood() {
        return 1 - subclonalLikelihood();
    }

}