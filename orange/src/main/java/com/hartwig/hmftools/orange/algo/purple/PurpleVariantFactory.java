package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;
import java.util.stream.Collectors;

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
import com.hartwig.hmftools.orange.algo.pave.PaveAlgo;
import com.hartwig.hmftools.orange.algo.pave.PaveEntry;
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.PurpleConversion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PurpleVariantFactory
{
    @NotNull
    private final PaveAlgo paveAlgo;

    public PurpleVariantFactory(@NotNull PaveAlgo paveAlgo)
    {
        this.paveAlgo = paveAlgo;
    }

    @Nullable
    public List<PurpleVariant> fromPurpleVariantContext(@Nullable List<PurpleVariantContext> contexts)
    {
        if(contexts == null)
        {
            return null;
        }
        return contexts.stream().map(this::fromPurpleVariantContext).collect(Collectors.toList());
    }

    @NotNull
    public PurpleVariant fromPurpleVariantContext(@NotNull PurpleVariantContext context)
    {
        List<PurpleTranscriptImpact> purpleVariantTranscriptImpacts =
                context.otherImpacts()
                        .stream()
                        .map(x -> PurpleConversion.convert(x, context.reportableTranscripts().contains(x.Transcript)))
                        .collect(Collectors.toList());
        PurpleAllelicDepth rnaDepth = context.rnaDepth() != null ? PurpleConversion.convert(context.rnaDepth()) : null;

        return ImmutablePurpleVariant.builder()
                .type(PurpleVariantType.valueOf(context.type().name()))
                .gene(context.gene())
                .chromosome(context.chromosome())
                .position(context.position())
                .ref(context.ref())
                .alt(context.alt())
                .worstCodingEffect(PurpleConversion.convert(context.worstCodingEffect()))
                .canonicalImpact(extractCanonicalImpact(context))
                .otherImpacts(purpleVariantTranscriptImpacts)
                .hotspot(HotspotType.valueOf(context.hotspot().name()))
                .tumorDepth(PurpleConversion.convert(context.allelicDepth()))
                .rnaDepth(rnaDepth)
                .adjustedCopyNumber(context.adjustedCopyNumber())
                .adjustedVAF(context.adjustedVAF())
                .minorAlleleCopyNumber(context.minorAlleleCopyNumber())
                .variantCopyNumber(context.variantCopyNumber())
                .biallelic(context.biallelic())
                .genotypeStatus(PurpleGenotypeStatus.valueOf(context.genotypeStatus().name()))
                .repeatCount(context.repeatCount())
                .subclonalLikelihood(context.subclonalLikelihood())
                .localPhaseSets(context.localPhaseSets())
                .build();
    }

    @NotNull
    private PurpleTranscriptImpact extractCanonicalImpact(PurpleVariantContext purpleContext)
    {
        PaveEntry paveEntry = paveAlgo.run(purpleContext.gene(), purpleContext.canonicalTranscript(), purpleContext.position());
        List<VariantEffect> variantEffects = VariantEffect.effectsToList(purpleContext.canonicalEffect());
        List<PurpleVariantEffect> purpleVariantEffects = ConversionUtil.mapToList(variantEffects, PurpleConversion::convert);
        return ImmutablePurpleTranscriptImpact.builder()
                .transcript(purpleContext.canonicalTranscript())
                .hgvsCodingImpact(purpleContext.canonicalHgvsCodingImpact())
                .hgvsProteinImpact(purpleContext.canonicalHgvsProteinImpact())
                .affectedCodon(paveEntry != null ? paveEntry.affectedCodon() : null)
                .affectedExon(paveEntry != null ? paveEntry.affectedExon() : null)
                .inSpliceRegion(purpleContext.spliceRegion())
                .effects(purpleVariantEffects)
                .codingEffect(PurpleConversion.convert(purpleContext.canonicalCodingEffect()))
                .reported(isCanonicalTranscriptReported(purpleContext))
                .build();
    }

    private boolean isCanonicalTranscriptReported(@NotNull PurpleVariantContext purpleContext)
    {
        if(purpleContext.reportableTranscripts().isEmpty())
        {
            return purpleContext.reported();
        }
        else
        {
            return purpleContext.reportableTranscripts().contains(purpleContext.canonicalTranscript());
        }
    }
}