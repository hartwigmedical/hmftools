package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.ImmutableAllelicDepthImpl;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.datamodel.purple.*;
import com.hartwig.hmftools.orange.algo.pave.PaveAlgo;
import com.hartwig.hmftools.orange.algo.pave.PaveEntry;
import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

public class PurpleVariantFactory {

    @NotNull
    private final PaveAlgo paveAlgo;

    public PurpleVariantFactory(@NotNull final PaveAlgo paveAlgo) {
        this.paveAlgo = paveAlgo;
    }

    @Nullable
    public List<PurpleVariant> create(@Nullable List<SomaticVariant> variants) {
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
    private PurpleVariant toPurpleVariant(@NotNull SomaticVariant variant) {
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
                .rnaDepth(Optional.ofNullable(variant.rnaDepth()).map(PurpleVariantFactory::convert).orElse(null))
                .referenceDepth(Optional.ofNullable(variant.referenceDepth()).map(PurpleVariantFactory::convert).orElse(null))
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
    private static PurpleAllelicDepth extractTumorDepth(@NotNull SomaticVariant variant) {
        return ImmutablePurpleAllelicDepth.builder()
                .alleleReadCount(variant.alleleReadCount())
                .totalReadCount(variant.totalReadCount())
                .build();
    }

    private static PurpleAllelicDepth convert(AllelicDepth allelicDepth) {
        return ImmutablePurpleAllelicDepth.builder()
                .alleleReadCount(allelicDepth.alleleReadCount())
                .totalReadCount(allelicDepth.totalReadCount())
                .build();
    }

    @NotNull
    private PurpleTranscriptImpact extractCanonicalImpact(@NotNull SomaticVariant variant) {
        // TODO Move effect parsing into SomaticVariant

        PaveEntry paveEntry = paveAlgo.run(variant.gene(), variant.canonicalTranscript(), variant.position());
        var variantEffects = VariantEffect.effectsToList(variant.canonicalEffect())
                .stream().map(effect -> PurpleVariantEffect.valueOf(effect.name()))
                .collect(Collectors.toList());
        return ImmutablePurpleTranscriptImpact.builder()
                .transcript(variant.canonicalTranscript())
                .hgvsCodingImpact(variant.canonicalHgvsCodingImpact())
                .hgvsProteinImpact(variant.canonicalHgvsProteinImpact())
                .affectedCodon(paveEntry != null ? paveEntry.affectedCodon() : null)
                .affectedExon(paveEntry != null ? paveEntry.affectedExon() : null)
                .spliceRegion(variant.spliceRegion())
                .effects(variantEffects)
                .codingEffect(PurpleCodingEffect.valueOf(variant.canonicalCodingEffect().name()))
                .build();
    }

    @NotNull
    private List<PurpleTranscriptImpact> extractOtherImpacts(@NotNull SomaticVariant variant) {
        List<PurpleTranscriptImpact> otherImpacts = Lists.newArrayList();
        // TODO Move other reported effects parsing into SomaticVariant
        // TODO Move effect parsing into SomaticVariant
        // TODO Add "splice region" details to non-canonical effects

        for (AltTranscriptReportableInfo altInfo : AltTranscriptReportableInfo.parseAltTranscriptInfo(variant.otherReportedEffects())) {
            PaveEntry paveEntry = paveAlgo.run(variant.gene(), altInfo.TransName, variant.position());
            var variantEffects = VariantEffect.effectsToList(altInfo.Effects)
                    .stream().map(effect -> PurpleVariantEffect.valueOf(effect.name()))
                    .collect(Collectors.toList());
            otherImpacts.add(ImmutablePurpleTranscriptImpact.builder()
                    .transcript(altInfo.TransName)
                    .hgvsCodingImpact(altInfo.HgvsCoding)
                    .hgvsProteinImpact(altInfo.HgvsProtein)
                    .affectedCodon(paveEntry != null ? paveEntry.affectedCodon() : null)
                    .affectedExon(paveEntry != null ? paveEntry.affectedExon() : null)
                    .spliceRegion(variant.spliceRegion())
                    .effects(variantEffects)
                    .codingEffect(PurpleCodingEffect.valueOf(altInfo.Effect.name()))
                    .build());
        }
        return otherImpacts;
    }
}
