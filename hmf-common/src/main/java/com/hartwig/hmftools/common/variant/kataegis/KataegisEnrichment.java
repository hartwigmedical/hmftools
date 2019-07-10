package com.hartwig.hmftools.common.variant.kataegis;

import java.util.function.Consumer;
import java.util.function.Predicate;

import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichment;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class KataegisEnrichment implements VariantContextEnrichment {

    public static final String KATAEGIS_FLAG = "KATAEGIS";
    private static final String KATAEGIS_FLAG_DESCRITION = "Strand kataegis detected on: FWD,REV or NONE";

    private final Predicate<VariantContext> forwardPredicate;
    private final Predicate<VariantContext> reversePredicate;

    private final KataegisQueue forwardDetector;
    private final KataegisQueue reverseDetector;

    public KataegisEnrichment(final Consumer<VariantContext> consumer) {

        forwardPredicate = x -> false;
        reversePredicate = x -> false;

        reverseDetector = new KataegisQueue(KataegisStatus.FWD, reversePredicate, consumer);
        forwardDetector = new KataegisQueue(KataegisStatus.FWD, forwardPredicate, reverseDetector);
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        forwardDetector.accept(context);
    }

    @Override
    public void flush() {
        forwardDetector.flush();
        reverseDetector.flush();
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template) {
        final VCFHeader outputVCFHeader = new VCFHeader(template.getMetaDataInInputOrder(), template.getSampleNamesInOrder());
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(KATAEGIS_FLAG,
                1,
                VCFHeaderLineType.String,
                KATAEGIS_FLAG_DESCRITION));

        return outputVCFHeader;
    }
}
