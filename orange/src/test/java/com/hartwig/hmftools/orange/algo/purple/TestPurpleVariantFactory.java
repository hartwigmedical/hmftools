package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.datamodel.purple.*;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class TestPurpleVariantFactory {

    private TestPurpleVariantFactory() {
    }

    @NotNull
    public static com.hartwig.hmftools.datamodel.purple.ImmutablePurpleVariant.Builder builder() {
        return com.hartwig.hmftools.datamodel.purple.ImmutablePurpleVariant.builder()
                .type(PurpleVariantType.UNDEFINED)
                .gene(Strings.EMPTY)
                .chromosome(Strings.EMPTY)
                .position(0)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .worstCodingEffect(PurpleCodingEffect.UNDEFINED)
                .canonicalImpact(impactBuilder().build())
                .hotspot(Hotspot.NON_HOTSPOT)
                .reported(false)
                .tumorDepth(depthBuilder().build())
                .rnaDepth(null)
                .adjustedCopyNumber(0)
                .adjustedVAF(0)
                .minorAlleleCopyNumber(0)
                .variantCopyNumber(0)
                .biallelic(false)
                .genotypeStatus(PurpleGenotypeStatus.UNKNOWN)
                .repeatCount(0)
                .subclonalLikelihood(0D)
                .localPhaseSets(null);
    }

    @NotNull
    public static ImmutablePurpleAllelicDepth.Builder depthBuilder() {
        return ImmutablePurpleAllelicDepth.builder().alleleReadCount(0).totalReadCount(0);
    }

    @NotNull
    public static com.hartwig.hmftools.datamodel.purple.ImmutablePurpleTranscriptImpact.Builder impactBuilder() {
        return com.hartwig.hmftools.datamodel.purple.ImmutablePurpleTranscriptImpact.builder()
                .transcript(Strings.EMPTY)
                .hgvsCodingImpact(Strings.EMPTY)
                .hgvsProteinImpact(Strings.EMPTY)
                .affectedCodon(null)
                .affectedExon(null)
                .spliceRegion(false)
                .codingEffect(PurpleCodingEffect.UNDEFINED);
    }
}
