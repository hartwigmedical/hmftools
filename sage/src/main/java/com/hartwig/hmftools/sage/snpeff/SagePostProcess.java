package com.hartwig.hmftools.sage.snpeff;

import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.CodingEffect.effect;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.MICROHOMOLOGY_FLAG;

import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.variant.CanonicalAnnotation;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory;
import com.hartwig.hmftools.sage.vcf.SageVCF;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

public class SagePostProcess implements AutoCloseable, Consumer<VariantContext> {

    private static final int BUFFER_DISTANCE = 10_000;
    static final String INFRAME_INSERTION = "inframe_insertion";
    static final String INFRAME_DELETION = "inframe_deletion";
    static final String SPLICE_DONOR_VARIANT = "splice_donor_variant";

    private final Consumer<VariantContext> consumer;
    private final List<VariantContext> buffer = Lists.newArrayList();
    private final CanonicalAnnotation canonicalAnnotation = new CanonicalAnnotation();
    private final Map<String, HmfTranscriptRegion> allGenesMap;

    public SagePostProcess(@NotNull final Map<String, HmfTranscriptRegion> allGenesMap, @NotNull final Consumer<VariantContext> consumer) {
        this.consumer = consumer;
        this.allGenesMap = allGenesMap;
    }

    @Override
    public void close() {
        flush();
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        if (!buffer.isEmpty()) {
            final VariantContext last = buffer.get(buffer.size() - 1);
            if (!last.getContig().equals(context.getContig()) || context.getStart() - last.getStart() > BUFFER_DISTANCE) {
                flush();
            }
        }

        process(context);
        buffer.add(context);
    }

    private void process(@NotNull final VariantContext context) {
        final List<SnpEffAnnotation> allAnnotations = SnpEffAnnotationFactory.fromContext(context);
        final Optional<SnpEffAnnotation> optionalCanonical = canonicalAnnotation.canonicalSnpEffAnnotation(allAnnotations);
        if (optionalCanonical.isPresent()) {
            final SnpEffAnnotation canonical = optionalCanonical.get();
            final CodingEffect codingEffect = effect(canonical.gene(), canonical.consequences());

            if (context.isSNP()) {
                spliceDonorPlus5(context, canonical);
            } else if (context.isIndel()) {
                phasedInframeIndel(context, canonical);
                inframeIndelWithMH(context, canonical, codingEffect);
            }
        }
    }

    private void phasedInframeIndel(@NotNull final VariantContext context, @NotNull final SnpEffAnnotation canonical) {
        int phase = context.getAttributeAsInt(SageVCF.PHASE, 0);
        if (phase == 0) {
            return;
        }

        final HmfTranscriptRegion transcript = allGenesMap.get(canonical.gene());
        if (transcript == null || context.getStart() < transcript.codingStart() || context.getEnd() > transcript.codingEnd()) {
            return;
        }

        final HmfExonRegion exon = selectExon(context, transcript);
        if (exon == null) {
            return;
        }

        int indelLength = indelLength(context);
        for (final VariantContext other : buffer) {
            if (other.isIndel() && other.getAttributeAsInt(SageVCF.PHASE, 0) == phase) {
                int otherLength = indelLength(other);
                if ((indelLength + otherLength) % 3 == 0 && other.getStart() >= exon.start() && other.getEnd() <= exon.end()) {
                    other.getCommonInfo().putAttribute(SagePostProcessVCF.PHASED_INFRAME_INDEL, true, true);
                    context.getCommonInfo().putAttribute(SagePostProcessVCF.PHASED_INFRAME_INDEL, true, true);
                }
            }
        }
    }

    @Nullable
    private HmfExonRegion selectExon(@NotNull final VariantContext context, @NotNull final HmfTranscriptRegion transcript) {
        for (HmfExonRegion exon : transcript.exome()) {
            if (context.getStart() >= exon.start() && context.getEnd() <= exon.end()) {
                return exon;
            }
        }

        return null;
    }

    private void spliceDonorPlus5(@NotNull final VariantContext context, @NotNull final SnpEffAnnotation canonical) {

        if (!context.isSNP() || !canonical.effects().contains("splice_region_variant")) {
            return;
        }

        final HmfTranscriptRegion transcript = allGenesMap.get(canonical.gene());
        if (transcript != null) {
            boolean forward = transcript.strand() == Strand.FORWARD;
            for (int i = 0; i < transcript.exome().size() - 1; i++) {

                final HmfExonRegion exon = transcript.exonByIndex(i);
                if (exon != null) {
                    if (forward && context.getStart() == exon.end() + 5 || !forward && context.getStart() == exon.start() - 5) {
                        insertEffectBeforeCanonicalAnnotation(context, canonical, SPLICE_DONOR_VARIANT);
                        return;
                    }
                }
            }
        }
    }

    private void inframeIndelWithMH(@NotNull final VariantContext context, @NotNull final SnpEffAnnotation canonical,
            @NotNull final CodingEffect codingEffect) {
        if (!context.isIndel() || !context.hasAttribute(MICROHOMOLOGY_FLAG) || !codingEffect.equals(NONSENSE_OR_FRAMESHIFT)
                || !canonical.effects().contains("splice") || !isMod3RefOrAlt(context)) {
            return;
        }

        int microhomologyLength = context.getAttributeAsString(MICROHOMOLOGY_FLAG, "").length();
        long rightAlignedIndelStart = context.getStart() + microhomologyLength + 1;
        long rightAlignedIndelEnd = context.getEnd() + microhomologyLength + 1;

        final HmfTranscriptRegion transcriptRegion = allGenesMap.get(canonical.gene());
        if (transcriptRegion != null && rightAlignedIndelStart >= transcriptRegion.codingStart() && rightAlignedIndelEnd <= transcriptRegion
                .codingEnd()) {
            for (int i = 1; i < transcriptRegion.exome().size(); i++) {

                final HmfExonRegion exon = transcriptRegion.exonByIndex(i);
                if (exon != null && rightAlignedIndelStart >= exon.start() && rightAlignedIndelEnd <= exon.end()) {

                    addInframeIndelEffectToCanonicalAnnotation(context, canonical);
                    return;
                }
            }
        }
    }

    private static void addInframeIndelEffectToCanonicalAnnotation(@NotNull final VariantContext context,
            @NotNull final SnpEffAnnotation canonical) {
        String effect = context.getReference().length() > context.getAlternateAllele(0).length() ? INFRAME_INSERTION : INFRAME_INSERTION;
        insertEffectBeforeCanonicalAnnotation(context, canonical, effect);
    }

    private static void insertEffectBeforeCanonicalAnnotation(@NotNull final VariantContext context,
            @NotNull final SnpEffAnnotation canonical, @NotNull final String newEffect) {

        final List<String> snpEffAnnotations = context.getAttributeAsStringList("ANN", "");
        for (int i = 0; i < snpEffAnnotations.size(); i++) {
            final String annotation = snpEffAnnotations.get(i);
            if (annotation.contains(canonical.effects()) && annotation.contains(canonical.transcript())) {
                final String updatedAnnotation = annotation.replace(canonical.effects(), newEffect);
                final String hmgTag = updatedAnnotation.charAt(updatedAnnotation.length() - 1) == '|' ? "HMF" : "&HMF";
                snpEffAnnotations.add(i, updatedAnnotation + hmgTag);
                context.getCommonInfo().putAttribute("ANN", snpEffAnnotations, true);
                return;
            }
        }
    }

    private int indelLength(@NotNull final VariantContext context) {
        if (context.isIndel()) {
            final String ref = context.getReference().getBaseString();
            final String alt = context.getAlternateAllele(0).getBaseString();
            return ref.length() > alt.length() ? -ref.length() + 1 : alt.length() - 1;
        }

        return 0;
    }

    private boolean isMod3RefOrAlt(@NotNull final VariantContext context) {
        int indelLength = Math.abs(indelLength(context));
        return (indelLength > 0 && indelLength % 3 == 0);
    }

    private void flush() {
        final Iterator<VariantContext> variantContextIterator = buffer.iterator();
        while (variantContextIterator.hasNext()) {
            consumer.accept(variantContextIterator.next());
            variantContextIterator.remove();
        }
    }
}
