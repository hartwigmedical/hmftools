package com.hartwig.hmftools.common.variant.enrich;

import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.variant.CanonicalAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class NeoEpitopeEnrichment implements VariantContextEnrichment {

    public static final String UP_STREAM_AMINO_ACIDS = "NEO_AAU";
    public static final String DOWN_STREAM_AMINO_ACIDS = "NEO_AAD";
    public static final String NOVEL_STREAM_AMINO_ACIDS = "NEO_AAN";
    public static final String NONSENSE_MEDIATED_DECAY = "NEO_NMD";

    private final IndexedFastaSequenceFile refGenome;
    private final CanonicalAnnotation canonicalAnnotationFactory;

    public NeoEpitopeEnrichment(final IndexedFastaSequenceFile refGenome) {
        this.refGenome = refGenome;
        this.canonicalAnnotationFactory = new CanonicalAnnotation();
    }

    @Override
    public void accept(@NotNull final VariantContext context) {

        final List<SnpEffAnnotation> allAnnotations = SnpEffAnnotationFactory.fromContext(context);
        final Optional<SnpEffAnnotation> canonicalAnnotation = canonicalAnnotationFactory.canonicalSnpEffAnnotation(allAnnotations);
        canonicalAnnotation.ifPresent(x -> enrich(x, context));
    }

    private void enrich(@NotNull final SnpEffAnnotation canonicalAnnotation, @NotNull final VariantContext context) {
        final String canonicalGene = canonicalAnnotation.gene();
        // TODO: do some neo-epitope logic
    }

    @Override
    public void flush() {
        // IGNORE
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template) {
        // TODO: Add the flags and descriptions here
        return template;
    }
}
