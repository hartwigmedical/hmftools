package com.hartwig.hmftools.orange.algo.purple;

import static com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo.parseAltTranscriptInfo;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.variant.SmallVariant;
import com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.datamodel.purple.HotspotType;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleVariant;
import com.hartwig.hmftools.datamodel.purple.PurpleAllelicDepth;
import com.hartwig.hmftools.datamodel.purple.PurpleGenotypeStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantType;
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.PurpleConversion;

import org.jetbrains.annotations.Nullable;

public final class PurpleVariantFactory
{
    @Nullable
    public static List<PurpleVariant> fromPurpleVariants(@Nullable final List<SmallVariant> variants, final List<DriverCatalog> drivers)
    {
        if(variants == null)
            return null;

        return variants.stream().map(x -> fromPurpleVariant(x, drivers)).collect(Collectors.toList());
    }

    public static PurpleVariant fromPurpleVariant(final SmallVariant variant, final List<DriverCatalog> drivers)
    {
        List<AltTranscriptReportableInfo> altTransEffects = parseAltTranscriptInfo(variant.otherReportedEffects());

        PurpleTranscriptImpact canonicalImpact = extractCanonicalImpact(variant);

        List<PurpleTranscriptImpact> nonCanonicalTransImpacts = Lists.newArrayList();

        for(AltTranscriptReportableInfo altTransInfo : altTransEffects)
        {
            // convert / filter etc
            if(variant.reportableTranscripts().contains(altTransInfo.TransName))
            {
                PurpleTranscriptImpact otherTransImpact = createOtherImpact(altTransInfo, canonicalImpact);
                nonCanonicalTransImpacts.add(otherTransImpact);
            }
        }

        PurpleAllelicDepth rnaDepth = variant.rnaDepth() != null ? PurpleConversion.convert(variant.rnaDepth()) : null;

        return ImmutablePurpleVariant.builder()
                .type(PurpleVariantType.valueOf(variant.type().name()))
                .gene(variant.gene())
                .chromosome(variant.chromosome())
                .position(variant.position())
                .ref(variant.ref())
                .alt(variant.alt())
                .worstCodingEffect(PurpleConversion.convert(variant.worstCodingEffect()))
                .canonicalImpact(canonicalImpact)
                .otherImpacts(nonCanonicalTransImpacts)
                .hotspot(HotspotType.valueOf(variant.hotspot().name()))
                .tumorDepth(PurpleConversion.convert(variant.allelicDepth()))
                .rnaDepth(rnaDepth)
                .adjustedCopyNumber(variant.adjustedCopyNumber())
                .adjustedVAF(variant.adjustedVAF())
                .minorAlleleCopyNumber(variant.minorAlleleCopyNumber())
                .variantCopyNumber(variant.variantCopyNumber())
                .biallelic(variant.biallelic())
                .biallelicProbability(variant.biallelicProbability())
                .genotypeStatus(PurpleGenotypeStatus.valueOf(variant.genotypeStatus().name()))
                .repeatCount(variant.repeatCount())
                .subclonalLikelihood(variant.subclonalLikelihood())
                .somaticLikelihood(PurpleConversion.convert(variant.somaticLikelihood()))
                .localPhaseSets(variant.localPhaseSets())
                .build();
    }

    private static PurpleTranscriptImpact createOtherImpact(
            final AltTranscriptReportableInfo transImpactInfo, final PurpleTranscriptImpact canonicalImpact)
    {
        List<VariantEffect> variantEffects = VariantEffect.effectsToList(transImpactInfo.Effects);
        List<PurpleVariantEffect> purpleVariantEffects = ConversionUtil.mapToList(variantEffects, PurpleConversion::convert);

        return ImmutablePurpleTranscriptImpact.builder()
                .transcript(transImpactInfo.TransName)
                .hgvsCodingImpact(transImpactInfo.HgvsCoding)
                .hgvsProteinImpact(transImpactInfo.HgvsProtein)
                .affectedCodon(canonicalImpact.affectedCodon())
                .affectedExon(canonicalImpact.affectedExon())
                .inSpliceRegion(canonicalImpact.inSpliceRegion())
                .effects(purpleVariantEffects)
                .codingEffect(canonicalImpact.codingEffect())
                .reported(true)
                .build();
    }

    private static PurpleTranscriptImpact extractCanonicalImpact(final SmallVariant variant)
    {
        List<VariantEffect> variantEffects = VariantEffect.effectsToList(variant.canonicalEffect());
        List<PurpleVariantEffect> purpleVariantEffects = ConversionUtil.mapToList(variantEffects, PurpleConversion::convert);

        return ImmutablePurpleTranscriptImpact.builder()
                .transcript(variant.canonicalTranscript())
                .hgvsCodingImpact(variant.canonicalHgvsCodingImpact())
                .hgvsProteinImpact(variant.canonicalHgvsProteinImpact())
                .affectedCodon(variant.canonicalAffectedCodon())
                .affectedExon(variant.canonicalAffectedExon())
                .inSpliceRegion(variant.spliceRegion())
                .effects(purpleVariantEffects)
                .codingEffect(PurpleConversion.convert(variant.canonicalCodingEffect()))
                .reported(variant.reported())
                .build();
    }
}