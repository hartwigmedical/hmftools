package com.hartwig.hmftools.orange.algo.purple;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.SUBCLONAL_LIKELIHOOD_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.filter.HumanChromosomeFilter;
import com.hartwig.hmftools.common.variant.filter.NTFilter;
import com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;
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
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.PurpleConversion;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;

public class PurpleVariantFactory {

    private final PaveAlgo paveAlgo;

    public PurpleVariantFactory(@NotNull PaveAlgo paveAlgo) {
        this.paveAlgo = paveAlgo;
    }

    @Nullable
    public List<PurpleVariant> fromPurpleVariantContext(@Nullable List<PurpleVariantContext> context) {
        if (context == null) {
            return null;
        }
        return context.stream().map(this::fromPurpleVariantContext).collect(Collectors.toList());
    }

    @NotNull
    public PurpleVariant fromPurpleVariantContext(@NotNull PurpleVariantContext context) {

        var purpleVariantTranscriptImpacts = context.otherImpacts().stream().map(PurpleConversion::convert).collect(Collectors.toList());
        var rnaDepth = context.rnaDepth() != null ? PurpleConversion.convert(context.rnaDepth()) : null;

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
                .hotspot(Hotspot.valueOf(context.hotspot().name()))
                .reported(context.reported())
                .tumorDepth(extractTumorDepth(context.tumorDepth()))
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

    private PurpleTranscriptImpact extractCanonicalImpact(PurpleVariantContext purpleContext) {
        var paveEntry = paveAlgo.run(purpleContext.gene(), purpleContext.canonicalTranscript(), purpleContext.position());
        List<VariantEffect> variantEffects = VariantEffect.effectsToList(purpleContext.canonicalEffect());
        List<PurpleVariantEffect> purpleVariantEffects = ConversionUtil.mapToList(variantEffects, PurpleConversion::convert);
        return ImmutablePurpleTranscriptImpact.builder()
                .transcript(purpleContext.canonicalTranscript())
                .hgvsCodingImpact(purpleContext.canonicalHgvsCodingImpact())
                .hgvsProteinImpact(purpleContext.canonicalHgvsProteinImpact())
                .affectedCodon(paveEntry != null ? paveEntry.affectedCodon() : null)
                .affectedExon(paveEntry != null ? paveEntry.affectedExon() : null)
                .spliceRegion(purpleContext.spliceRegion())
                .effects(purpleVariantEffects)
                .codingEffect(PurpleConversion.convert(purpleContext.canonicalCodingEffect()))
                .build();
    }

    private static PurpleAllelicDepth extractTumorDepth(AllelicDepth tumorDepth) {
        return ImmutablePurpleAllelicDepth.builder()
                .alleleReadCount(tumorDepth.alleleReadCount())
                .totalReadCount(tumorDepth.totalReadCount())
                .build();
    }
}