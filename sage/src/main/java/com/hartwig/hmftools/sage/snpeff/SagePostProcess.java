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

import org.jetbrains.annotations.NotNull;

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

            spliceDonorPlus5(context, canonical);
            inframeIndelWithMH(context, canonical, codingEffect);

        }
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

    private static void insertEffectBeforeCanonicalAnnotation(@NotNull final VariantContext context, @NotNull final SnpEffAnnotation canonical,
            @NotNull final String newEffect) {

        final List<String> snpEffAnnotations = context.getAttributeAsStringList("ANN", "");
        for (int i = 0; i < snpEffAnnotations.size(); i++) {
            final String annotation = snpEffAnnotations.get(i);
            if (annotation.contains(canonical.effects()) && annotation.contains(canonical.transcript())) {
                final String updatedAnnotation = annotation.replace(canonical.effects(), newEffect);
                final String hmgTag = updatedAnnotation.charAt(updatedAnnotation.length() -1) == '|' ? "HMF" : "&HMF";
                snpEffAnnotations.add(i, updatedAnnotation + hmgTag);
                context.getCommonInfo().putAttribute("ANN", snpEffAnnotations, true);
                return;
            }
        }
    }

    private boolean isMod3RefOrAlt(@NotNull final VariantContext context) {
        int ref = context.getReference().getBaseString().length() - 1;
        int alt = context.getReference().getBaseString().length() - 1;

        return (ref > 0 && ref % 3 == 0) || (alt > 0 && alt % 3 == 0);
    }

    private void flush() {
        final Iterator<VariantContext> variantContextIterator = buffer.iterator();
        while (variantContextIterator.hasNext()) {
            consumer.accept(variantContextIterator.next());
            variantContextIterator.remove();
        }
    }
}
