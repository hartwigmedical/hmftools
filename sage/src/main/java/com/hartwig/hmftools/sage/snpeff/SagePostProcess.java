package com.hartwig.hmftools.sage.snpeff;

import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.CodingEffect.effect;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.MICROHOMOLOGY_FLAG;

import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.CanonicalTranscript;
import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.sage.SagePostProcessVCF;
import com.hartwig.hmftools.common.variant.CanonicalAnnotation;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantConsequence;
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

    private final List<VariantContext> buffer;
    private final Consumer<VariantContext> consumer;
    private final CanonicalAnnotation canonicalAnnotation;
    private final Map<String, HmfTranscriptRegion> allGenesMap;

    public SagePostProcess(@NotNull final Map<String, HmfTranscriptRegion> allGenesMap,
            @NotNull final List<CanonicalTranscript> transcripts, @NotNull final Consumer<VariantContext> consumer) {
        this.buffer = Lists.newArrayList();
        this.consumer = consumer;
        this.allGenesMap = allGenesMap;
        this.canonicalAnnotation = new CanonicalAnnotation(transcripts);
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

            context.getCommonInfo().putAttribute(SagePostProcessVCF.HMF_CANONICAL_GENE, canonical.gene());
            context.getCommonInfo().putAttribute(SagePostProcessVCF.HMF_CANONICAL_EFFECT, consequences(canonical));
            context.getCommonInfo().putAttribute(SagePostProcessVCF.HMF_CANONICAL_IMPACT, codingEffect.toString());

            if (spliceDonorPlus5(context, canonical)) {
                context.getCommonInfo().putAttribute(SagePostProcessVCF.HMF_CANONICAL_EFFECT, "splice_donor_5_base_variant", true);
                context.getCommonInfo().putAttribute(SagePostProcessVCF.HMF_CANONICAL_IMPACT, CodingEffect.SPLICE.toString(), true);
            } else if (inframeIndelWithMH(context, canonical, codingEffect)) {
                context.getCommonInfo().putAttribute(SagePostProcessVCF.HMF_CANONICAL_EFFECT, "inframe_mh", true);
                context.getCommonInfo().putAttribute(SagePostProcessVCF.HMF_CANONICAL_IMPACT, CodingEffect.MISSENSE.toString(), true);
            }

            processPhasedInframeIndel(context, canonical);
            processMixedGermlineImpact(context);
        }
    }

    @NotNull
    private List<String> consequences(@NotNull final SnpEffAnnotation annotation) {
        return annotation.consequences().stream().map(VariantConsequence::parentSequenceOntologyTerm).collect(Collectors.toList());
    }

    private void processMixedGermlineImpact(@NotNull final VariantContext context) {
        int mixedImpact = context.getAttributeAsInt(SageVCF.MIXED_GERMLINE_IMPACT, 0);
        if (mixedImpact == 0) {
            return;
        }

        boolean newIsSnv = context.isSNP();
        boolean newIsMnv = context.isMNP();

        for (VariantContext otherContext : buffer) {
            if (otherContext.getAttributeAsInt(SageVCF.MIXED_GERMLINE_IMPACT, 0) == mixedImpact) {

                boolean otherIsSnv = otherContext.isSNP();
                boolean otherIsMnv = otherContext.isMNP();

                final VariantContext snv;
                final VariantContext mnv;
                if (newIsSnv && otherIsMnv) {
                    snv = context;
                    mnv = otherContext;
                } else if (newIsMnv && otherIsSnv) {
                    snv = otherContext;
                    mnv = context;
                } else {
                    snv = null;
                    mnv = null;
                }

                if (snv != null && mnv != null) {
                    final List<SnpEffAnnotation> mnvAnnotations = SnpEffAnnotationFactory.fromContext(context);
                    final Optional<SnpEffAnnotation> optionalCanonical = canonicalAnnotation.canonicalSnpEffAnnotation(mnvAnnotations);
                    if (optionalCanonical.isPresent()) {
                        final SnpEffAnnotation canonical = optionalCanonical.get();
                        final CodingEffect codingEffect = effect(canonical.gene(), canonical.consequences());

                        snv.getCommonInfo().putAttribute(SagePostProcessVCF.HMF_CANONICAL_GENE, canonical.gene(), true);
                        snv.getCommonInfo().putAttribute(SagePostProcessVCF.HMF_CANONICAL_EFFECT, consequences(canonical), true);
                        snv.getCommonInfo().putAttribute(SagePostProcessVCF.HMF_CANONICAL_IMPACT, codingEffect.toString(), true);
                    }
                }
            }
        }
    }

    private void processPhasedInframeIndel(@NotNull final VariantContext context, @NotNull final SnpEffAnnotation canonical) {
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

                    other.getCommonInfo().putAttribute(SagePostProcessVCF.HMF_CANONICAL_EFFECT, "inframe_phased", true);
                    context.getCommonInfo().putAttribute(SagePostProcessVCF.HMF_CANONICAL_EFFECT, "inframe_phased", true);

                    other.getCommonInfo().putAttribute(SagePostProcessVCF.HMF_CANONICAL_IMPACT, CodingEffect.MISSENSE.toString(), true);
                    context.getCommonInfo().putAttribute(SagePostProcessVCF.HMF_CANONICAL_IMPACT, CodingEffect.MISSENSE.toString(), true);
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

    private boolean spliceDonorPlus5(@NotNull final VariantContext context, @NotNull final SnpEffAnnotation canonical) {

        if (!context.isSNP() || !canonical.effects().contains("splice_region_variant")) {
            return false;
        }

        final HmfTranscriptRegion transcript = allGenesMap.get(canonical.gene());
        if (transcript != null) {
            boolean forward = transcript.strand() == Strand.FORWARD;
            for (int i = 0; i < transcript.exome().size() - 1; i++) {

                final HmfExonRegion exon = transcript.exonByIndex(i);
                if (exon != null) {
                    if (forward && context.getStart() == exon.end() + 5 || !forward && context.getStart() == exon.start() - 5) {
                        return true;
                    }
                }
            }
        }

        return false;
    }

    private boolean inframeIndelWithMH(@NotNull final VariantContext context, @NotNull final SnpEffAnnotation canonical,
            @NotNull final CodingEffect codingEffect) {
        if (!context.isIndel() || !context.hasAttribute(MICROHOMOLOGY_FLAG) || !codingEffect.equals(NONSENSE_OR_FRAMESHIFT)
                || !canonical.effects().contains("splice") || !isMod3RefOrAlt(context)) {
            return false;
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
                    return true;
                }
            }
        }

        return false;
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
