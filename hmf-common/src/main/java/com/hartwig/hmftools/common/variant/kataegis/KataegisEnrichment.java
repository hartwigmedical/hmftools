package com.hartwig.hmftools.common.variant.kataegis;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.variant.RefContextEnrichment;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichment;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class KataegisEnrichment implements VariantContextEnrichment {

    public static final String KATAEGIS_FLAG = "KATAEGIS";
    private static final String KATAEGIS_FLAG_DESCRITION = "Strand kataegis detected on: FWD,REV or NONE";

    private final KataegisQueue forwardDetector;
    private final KataegisQueue reverseDetector;

    public KataegisEnrichment(final Consumer<VariantContext> consumer) {
        reverseDetector = new KataegisQueue(KataegisStatus.REV, KataegisEnrichment::isReverseCandidate, consumer);
        forwardDetector = new KataegisQueue(KataegisStatus.FWD, KataegisEnrichment::isForwardCandidate, reverseDetector);
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
        template.addMetaDataLine(new VCFInfoHeaderLine(KATAEGIS_FLAG, 1, VCFHeaderLineType.String, KATAEGIS_FLAG_DESCRITION));

        return template;
    }

    private static boolean isForwardCandidate(@NotNull final VariantContext context) {

        final boolean altMatch =
                context.getAlternateAlleles().stream().anyMatch(x -> x.getBaseString().equals("T") || x.getBaseString().equals("G"));

        final String triContext = context.getAttributeAsString(RefContextEnrichment.TRINUCLEOTIDE_FLAG, Strings.EMPTY);
        final boolean triMatch = triContext.startsWith("TC");

        return triMatch && altMatch;
    }

    private static boolean isReverseCandidate(@NotNull final VariantContext context) {

        final boolean altMatch =
                context.getAlternateAlleles().stream().anyMatch(x -> x.getBaseString().equals("C") || x.getBaseString().equals("A"));

        final String triContext = context.getAttributeAsString(RefContextEnrichment.TRINUCLEOTIDE_FLAG, Strings.EMPTY);
        final boolean triMatch = triContext.endsWith("GA");

        return triMatch && altMatch;

    }

}
