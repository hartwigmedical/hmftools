package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.datamodel.purple.HotspotType;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleAllelicDepth;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleVariant;
import com.hartwig.hmftools.datamodel.purple.PurpleCodingEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleGenotypeStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantType;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class TestPurpleVariantFactory
{
    @NotNull
    public static ImmutablePurpleVariant.Builder builder()
    {
        return ImmutablePurpleVariant.builder()
                .type(PurpleVariantType.UNDEFINED)
                .gene(Strings.EMPTY)
                .chromosome(Strings.EMPTY)
                .position(0)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .worstCodingEffect(PurpleCodingEffect.UNDEFINED)
                .canonicalImpact(impactBuilder().build())
                .hotspot(HotspotType.NON_HOTSPOT)
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
    public static ImmutablePurpleAllelicDepth.Builder depthBuilder()
    {
        return ImmutablePurpleAllelicDepth.builder().alleleReadCount(0).totalReadCount(0);
    }

    @NotNull
    public static ImmutablePurpleTranscriptImpact.Builder impactBuilder()
    {
        return ImmutablePurpleTranscriptImpact.builder()
                .transcript(Strings.EMPTY)
                .hgvsCodingImpact(Strings.EMPTY)
                .hgvsProteinImpact(Strings.EMPTY)
                .affectedCodon(null)
                .affectedExon(null)
                .inSpliceRegion(false)
                .codingEffect(PurpleCodingEffect.UNDEFINED);
    }
}
