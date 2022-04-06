package com.hartwig.hmftools.common.test;

import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class SomaticVariantTestBuilderFactory {

    private SomaticVariantTestBuilderFactory() {
    }

    @NotNull
    public static ImmutableSomaticVariantImpl.Builder create() {
        return ImmutableSomaticVariantImpl.builder()
                .qual(100)
                .chromosome(Strings.EMPTY)
                .position(0)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .type(VariantType.UNDEFINED)
                .filter(Strings.EMPTY)
                .totalReadCount(0)
                .alleleReadCount(0)
                .gene(Strings.EMPTY)
                .genotypeStatus(GenotypeStatus.UNKNOWN)
                .genesAffected(0)
                .canonicalEffect(Strings.EMPTY)
                .canonicalTranscript(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.UNDEFINED)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .spliceRegion(false)
                .otherReportedEffects(Strings.EMPTY)
                .worstCodingEffect(CodingEffect.NONE)
                .hotspot(Hotspot.NON_HOTSPOT)
                .recovered(false)
                .reported(false)
                .adjustedCopyNumber(0D)
                .adjustedVAF(0D)
                .minorAlleleCopyNumber(0D)
                .germlineStatus(GermlineStatus.UNKNOWN)
                .variantCopyNumber(0)
                .biallelic(false)
                .kataegis(Strings.EMPTY)
                .trinucleotideContext(Strings.EMPTY)
                .microhomology(Strings.EMPTY)
                .repeatSequence(Strings.EMPTY)
                .repeatCount(0)
                .subclonalLikelihood(0)
                .tier(VariantTier.UNKNOWN)
                .mappability(0D);
    }
}