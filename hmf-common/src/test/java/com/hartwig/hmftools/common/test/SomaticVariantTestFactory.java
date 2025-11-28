package com.hartwig.hmftools.common.test;

import static com.hartwig.hmftools.common.variant.SomaticLikelihood.UNKNOWN;

import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;
import com.hartwig.hmftools.common.variant.ImmutableVariantImpl;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class SomaticVariantTestFactory
{
    @NotNull
    public static ImmutableSomaticVariantImpl.Builder builder()
    {
        return ImmutableSomaticVariantImpl.builder()
                .variant(variantBuilder().build())
                .kataegis(Strings.EMPTY)
                .subclonalLikelihood(0)
                .gnomadFrequency(0)
                .somaticLikelihood(UNKNOWN);
    }

    @NotNull static ImmutableVariantImpl.Builder variantBuilder()
    {
        return ImmutableVariantImpl.builder()
                .chromosome(Strings.EMPTY)
                .position(0)
                .allelicDepth(new AllelicDepth(0, 0))
                .type(VariantType.UNDEFINED)
                .gene(Strings.EMPTY)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .canonicalTranscript(Strings.EMPTY)
                .canonicalEffect(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.UNDEFINED)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .qual(100)
                .mappability(0)
                .filter(Strings.EMPTY)
                .genesAffected(0)
                .spliceRegion(false)
                .otherReportedEffects(Strings.EMPTY)
                .worstCodingEffect(CodingEffect.UNDEFINED)
                .tier(VariantTier.UNKNOWN)
                .hotspot(Hotspot.NON_HOTSPOT)
                .reported(false)
                .adjustedCopyNumber(0)
                .adjustedVAF(0)
                .minorAlleleCopyNumber(0)
                .variantCopyNumber(0)
                .biallelic(false)
                .genotypeStatus(GenotypeStatus.UNKNOWN)
                .germlineStatus(GermlineStatus.UNKNOWN)
                .trinucleotideContext(Strings.EMPTY)
                .microhomology(Strings.EMPTY)
                .repeatSequence(Strings.EMPTY)
                .repeatCount(0);
    }
}
