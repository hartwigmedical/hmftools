package com.hartwig.hmftools.common.variant.enrich;

import java.util.Set;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.variant.kataegis.KataegisQueue;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class KataegisEnrichment implements VariantContextEnrichment {

    public static final String KATAEGIS_FLAG = "KT";
    private static final String KATAEGIS_FLAG_DESCRIPTION = "Forward/reverse kataegis id";

    private final KataegisQueue forwardDetector;
    private final KataegisQueue reverseDetector;

    public KataegisEnrichment(@NotNull final Consumer<VariantContext> consumer) {
        reverseDetector = new KataegisQueue("REV", KataegisEnrichment::isReverseCandidate, consumer);
        forwardDetector = new KataegisQueue("FWD", KataegisEnrichment::isForwardCandidate, reverseDetector::accept);
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
        template.addMetaDataLine(new VCFInfoHeaderLine(KATAEGIS_FLAG, 1, VCFHeaderLineType.String, KATAEGIS_FLAG_DESCRIPTION));

        return template;
    }

    private static boolean isForwardCandidate(@NotNull final VariantContext context) {
        final boolean altMatch =
                context.getAlternateAlleles().stream().anyMatch(x -> x.getBaseString().equals("T") || x.getBaseString().equals("G"));

        final String triContext = context.getAttributeAsString(SomaticRefContextEnrichment.TRINUCLEOTIDE_FLAG, Strings.EMPTY);
        final boolean triMatch = triContext.startsWith("TC");

        return isNotFiltered(context) && triMatch && altMatch;
    }

    private static boolean isReverseCandidate(@NotNull final VariantContext context) {
        final boolean altMatch =
                context.getAlternateAlleles().stream().anyMatch(x -> x.getBaseString().equals("C") || x.getBaseString().equals("A"));

        final String triContext = context.getAttributeAsString(SomaticRefContextEnrichment.TRINUCLEOTIDE_FLAG, Strings.EMPTY);
        final boolean triMatch = triContext.endsWith("GA");

        return isNotFiltered(context) && triMatch && altMatch;
    }

    private static boolean isNotFiltered(@NotNull final VariantContext context) {
        final Set<String> filters = context.getFilters();
        return filters.isEmpty() || (filters.size() == 1 && filters.contains("PASS"));
    }
}
