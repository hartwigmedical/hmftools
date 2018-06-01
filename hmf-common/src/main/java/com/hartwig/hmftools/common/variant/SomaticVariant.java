package com.hartwig.hmftools.common.variant;

import java.util.List;

import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.variant.snpeff.VariantAnnotation;

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

    @Nullable
    String cosmicID();

    @NotNull
    List<VariantAnnotation> annotations();

    @NotNull
    String gene();

    int genesEffected();

    @NotNull
    String worstEffect();

    @NotNull
    String worstEffectTranscript();

    @NotNull
    CodingEffect worstCodingEffect();

    boolean hotspot();

    double mappability();

    default boolean isDBSNP() {
        return dbsnpID() != null;
    }

    default boolean isCOSMIC() {
        return cosmicID() != null;
    }

    default boolean isFiltered() {return !filter().equals("PASS");}
}