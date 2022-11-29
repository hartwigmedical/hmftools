package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableAllelicDepthImpl;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class PurpleVariantTestFactory {

    private PurpleVariantTestFactory() {
    }

    @NotNull
    public static ImmutablePurpleVariant.Builder builder() {
        return ImmutablePurpleVariant.builder()
                .type(VariantType.UNDEFINED)
                .gene(Strings.EMPTY)
                .genesAffected(0)
                .chromosome(Strings.EMPTY)
                .position(0)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .worstCodingEffect(CodingEffect.UNDEFINED)
                .canonicalImpact(impactBuilder().build())
                .hotspot(Hotspot.NON_HOTSPOT)
                .reported(false)
                .filtered(false)
                .filter(Strings.EMPTY)
                .recovered(false)
                .tumorDepth(depthBuilder().build())
                .rnaDepth(null)
                .referenceDepth(null)
                .adjustedCopyNumber(0)
                .adjustedVAF(0)
                .minorAlleleCopyNumber(0)
                .variantCopyNumber(0)
                .biallelic(false)
                .genotypeStatus(GenotypeStatus.UNKNOWN)
                .germlineStatus(GermlineStatus.UNKNOWN)
                .trinucleotideContext(Strings.EMPTY)
                .mappability(0)
                .microhomology(Strings.EMPTY)
                .repeatSequence(Strings.EMPTY)
                .repeatCount(0)
                .kataegis(Strings.EMPTY)
                .tier(VariantTier.UNKNOWN)
                .subclonalLikelihood(0D)
                .localPhaseSets(null);
    }

    @NotNull
    public static ImmutableAllelicDepthImpl.Builder depthBuilder() {
        return ImmutableAllelicDepthImpl.builder().alleleReadCount(0).totalReadCount(0);
    }

    @NotNull
    public static ImmutablePurpleTranscriptImpact.Builder impactBuilder() {
        return ImmutablePurpleTranscriptImpact.builder()
                .transcriptId(Strings.EMPTY)
                .hgvsCodingImpact(Strings.EMPTY)
                .hgvsProteinImpact(Strings.EMPTY)
                .affectedCodon(null)
                .affectedExon(null)
                .spliceRegion(false)
                .codingEffect(CodingEffect.UNDEFINED);
    }
}
