package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;

import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.datamodel.purple.Hotspot;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleAllelicDepth;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleVariant;
import com.hartwig.hmftools.datamodel.purple.PurpleAllelicDepth;
import com.hartwig.hmftools.datamodel.purple.PurpleGenotypeStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantType;
import com.hartwig.hmftools.orange.algo.pave.PaveAlgo;
import com.hartwig.hmftools.orange.algo.pave.PaveEntry;
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.PurpleConversion;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

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
        com.hartwig.hmftools.common.variant.AllelicDepth nullable = variant.rnaDepth();
        return ImmutablePurpleVariant.builder()
                .type(PurpleVariantType.valueOf(variant.type().name()))
                .gene(variant.gene())
                .chromosome(variant.chromosome())
                .position(variant.position())
                .ref(variant.ref())
                .alt(variant.alt())
                .worstCodingEffect(PurpleConversion.convert(variant.worstCodingEffect()))
                .canonicalImpact(extractCanonicalImpact(variant))
                .otherImpacts(extractOtherImpacts(variant))
                .hotspot(Hotspot.valueOf(variant.hotspot().name()))
                .reported(variant.reported())
                .tumorDepth(extractTumorDepth(variant))
                .rnaDepth(nullable == null ? null : PurpleConversion.convert(nullable))
                .adjustedCopyNumber(variant.adjustedCopyNumber())
                .adjustedVAF(variant.adjustedVAF())
                .minorAlleleCopyNumber(variant.minorAlleleCopyNumber())
                .variantCopyNumber(variant.variantCopyNumber())
                .biallelic(variant.biallelic())
                .genotypeStatus(PurpleGenotypeStatus.valueOf(variant.genotypeStatus().name()))
                .repeatCount(variant.repeatCount())
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

    @NotNull
    private PurpleTranscriptImpact extractCanonicalImpact(@NotNull SomaticVariant variant)
    {
        // TODO Move effect parsing into SomaticVariant

        PaveEntry paveEntry = paveAlgo.run(variant.gene(), variant.canonicalTranscript(), variant.position());
        List<VariantEffect> variantEffects = VariantEffect.effectsToList(variant.canonicalEffect());
        List<PurpleVariantEffect> purpleVariantEffects = ConversionUtil.mapCollection(variantEffects, PurpleConversion::convert);
        return ImmutablePurpleTranscriptImpact.builder()
                .transcript(variant.canonicalTranscript())
                .hgvsCodingImpact(variant.canonicalHgvsCodingImpact())
                .hgvsProteinImpact(variant.canonicalHgvsProteinImpact())
                .affectedCodon(paveEntry != null ? paveEntry.affectedCodon() : null)
                .affectedExon(paveEntry != null ? paveEntry.affectedExon() : null)
                .spliceRegion(variant.spliceRegion())
                .effects(purpleVariantEffects)
                .codingEffect(PurpleConversion.convert(variant.canonicalCodingEffect()))
                .build();
    }

    @NotNull
    private List<PurpleTranscriptImpact> extractOtherImpacts(@NotNull SomaticVariant variant) {
        List<PurpleTranscriptImpact> otherImpacts = Lists.newArrayList();
        // TODO Move other reported effects parsing into SomaticVariant
        // TODO Move effect parsing into SomaticVariant
        // TODO Add "splice region" details to non-canonical effects

        for (AltTranscriptReportableInfo altInfo : AltTranscriptReportableInfo.parseAltTranscriptInfo(variant.otherReportedEffects()))
        {
            PaveEntry paveEntry = paveAlgo.run(variant.gene(), altInfo.TransName, variant.position());
            List<VariantEffect> variantEffects = VariantEffect.effectsToList(altInfo.Effects);
            List<PurpleVariantEffect> purpleVariantEffects = ConversionUtil.mapCollection(variantEffects, PurpleConversion::convert);
            otherImpacts.add(ImmutablePurpleTranscriptImpact.builder()
                    .transcript(altInfo.TransName)
                    .hgvsCodingImpact(altInfo.HgvsCoding)
                    .hgvsProteinImpact(altInfo.HgvsProtein)
                    .affectedCodon(paveEntry != null ? paveEntry.affectedCodon() : null)
                    .affectedExon(paveEntry != null ? paveEntry.affectedExon() : null)
                    .spliceRegion(variant.spliceRegion())
                    .effects(purpleVariantEffects)
                    .codingEffect(PurpleConversion.convert(altInfo.Effect))
                    .build());
        }
        return otherImpacts;
    }
}
