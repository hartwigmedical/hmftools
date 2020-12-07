package com.hartwig.hmftools.common.pathogenic;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class PathogenicSummaryFactory {

    public static final String CLNSIG = "CLNSIG";
    public static final String CLNSIGCONF = "CLNSIGCONF";
    public static final String BLACKLIST_BED = "BLACKLIST_BED";
    public static final String BLACKLIST_VCF = "BLACKLIST_VCF";

    @NotNull
    public static PathogenicSummary fromContext(@NotNull VariantContext context) {

        final String clnSig = clnSig(context);
        final String clnSigConf = clnSigConf(context);
        final String clinvarInfo = !clnSigConf.isEmpty() ? clnSigConf : clnSig;

        final Pathogenic path = context.getAttributeAsBoolean(BLACKLIST_BED, false) || context.getAttributeAsBoolean(BLACKLIST_VCF, false)
                ? Pathogenic.BENIGN_BLACKLIST
                : Pathogenic.fromClinvarAnnotation(clnSig, clnSigConf);

        return ImmutablePathogenicSummary.builder().clinvarInfo(clinvarInfo).pathogenicity(path).build();
    }

    public static String clnSig(@NotNull VariantContext context) {
        return context.getAttributeAsString(CLNSIG, Strings.EMPTY);
    }

    public static String clnSigConf(@NotNull VariantContext context) {
        return String.join(",", context.getAttributeAsStringList(CLNSIGCONF, Strings.EMPTY));
    }

}
