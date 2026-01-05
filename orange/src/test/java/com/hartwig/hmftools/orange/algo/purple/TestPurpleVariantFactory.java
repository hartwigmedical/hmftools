package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableSmallVariantImpl;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
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
                .tumorDepth(depthBuilder().build())
                .rnaDepth(null)
                .adjustedCopyNumber(0)
                .adjustedVAF(0)
                .minorAlleleCopyNumber(0)
                .variantCopyNumber(0)
                .biallelic(false)
                .biallelicProbability(0.)
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
                .codingEffect(PurpleCodingEffect.UNDEFINED)
                .reported(false);
    }

    public static ImmutablePurpleVariantContext.Builder contextBuilder()
    {
        return ImmutablePurpleVariantContext.builder()
                .from(variantBuilder().build())
                .biallelicProbability(0)
                .subclonalLikelihood(0);
    }

    @NotNull static ImmutableSmallVariantImpl.Builder variantBuilder()
    {
        return ImmutableSmallVariantImpl.builder()
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
                .qual(0)
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
                .repeatCount(0)
                .subclonalLikelihood(0)
                .gnomadFrequency(0)
                .pathogenic(false);
    }
}
