package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;

import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.ImmutableAllelicDepthImpl;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class PurpleVariantFactory {

    private PurpleVariantFactory() {
    }

    @Nullable
    public static List<PurpleVariant> create(@Nullable List<SomaticVariant> variants) {
        if (variants == null) {
            return null;
        }

        List<PurpleVariant> purpleVariants = Lists.newArrayList();
        for (SomaticVariant variant : variants) {
            purpleVariants.add(toPurpleVariant(variant));
        }
        return purpleVariants;
    }

    @NotNull
    private static PurpleVariant toPurpleVariant(@NotNull SomaticVariant variant) {
        return ImmutablePurpleVariant.builder()
                .type(variant.type())
                .gene(variant.gene())
                .genesAffected(variant.genesAffected())
                .chromosome(variant.chromosome())
                .position(variant.position())
                .ref(variant.ref())
                .alt(variant.alt())
                .worstCodingEffect(variant.worstCodingEffect())
                .canonicalImpact(extractCanonicalImpact(variant))
                .otherImpacts(extractOtherImpacts(variant))
                .hotspot(variant.hotspot())
                .reported(variant.reported())
                .filtered(variant.isFiltered())
                .filter(variant.filter())
                .recovered(variant.recovered())
                .tumorDepth(extractTumorDepth(variant))
                .rnaDepth(variant.rnaDepth())
                .referenceDepth(variant.referenceDepth())
                .adjustedCopyNumber(variant.adjustedCopyNumber())
                .adjustedVAF(variant.adjustedVAF())
                .minorAlleleCopyNumber(variant.minorAlleleCopyNumber())
                .variantCopyNumber(variant.variantCopyNumber())
                .biallelic(variant.biallelic())
                .genotypeStatus(variant.genotypeStatus())
                .germlineStatus(variant.germlineStatus())
                .trinucleotideContext(variant.trinucleotideContext())
                .mappability(variant.mappability())
                .microhomology(variant.microhomology())
                .repeatSequence(variant.repeatSequence())
                .repeatCount(variant.repeatCount())
                .kataegis(variant.kataegis())
                .tier(variant.tier())
                .subclonalLikelihood(variant.subclonalLikelihood())
                .localPhaseSets(variant.localPhaseSets())
                .build();
    }

    @NotNull
    private static AllelicDepth extractTumorDepth(@NotNull SomaticVariant variant) {
        return ImmutableAllelicDepthImpl.builder()
                .alleleReadCount(variant.alleleReadCount())
                .totalReadCount(variant.totalReadCount())
                .build();
    }

    @NotNull
    private static PurpleTranscriptImpact extractCanonicalImpact(@NotNull SomaticVariant variant) {
        // TODO Populate codon and exon.
        // TODO Move effect parsing into SomaticVariant

        return ImmutablePurpleTranscriptImpact.builder()
                .transcript(variant.canonicalTranscript())
                .hgvsCodingImpact(variant.canonicalHgvsCodingImpact())
                .hgvsProteinImpact(variant.canonicalHgvsProteinImpact())
                .affectedCodon(null)
                .affectedExon(null)
                .spliceRegion(variant.spliceRegion())
                .effects(VariantEffect.effectsToList(variant.canonicalEffect()))
                .codingEffect(variant.canonicalCodingEffect())
                .build();
    }

    @NotNull
    private static List<PurpleTranscriptImpact> extractOtherImpacts(@NotNull SomaticVariant variant) {
        List<PurpleTranscriptImpact> otherImpacts = Lists.newArrayList();
        for (AltTranscriptReportableInfo altInfo : AltTranscriptReportableInfo.parseAltTranscriptInfo(variant.otherReportedEffects())) {
            otherImpacts.add(ImmutablePurpleTranscriptImpact.builder()
                    .transcript(altInfo.TransName)
                    .hgvsCodingImpact(altInfo.HgvsCoding)
                    .hgvsProteinImpact(altInfo.HgvsProtein)
                    .affectedCodon(null)
                    .affectedExon(null)
                    .spliceRegion(variant.spliceRegion())
                    .effects(VariantEffect.effectsToList(altInfo.Effects))
                    .codingEffect(altInfo.Effect)
                    .build());
        }
        return otherImpacts;
    }
}
