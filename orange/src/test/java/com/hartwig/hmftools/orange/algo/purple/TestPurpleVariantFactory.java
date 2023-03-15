package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.variant.*;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleAllelicDepth;
import com.hartwig.hmftools.datamodel.purple.PurpleCodingEffect;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class TestPurpleVariantFactory {

    private TestPurpleVariantFactory() {
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
