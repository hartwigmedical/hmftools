package com.hartwig.hmftools.common.pathogenic;

import org.apache.logging.log4j.util.Strings;

import htsjdk.variant.variantcontext.VariantContext;

public final class PathogenicSummaryFactory
{
    public static final String CLNSIG = "CLNSIG";
    public static final String CLNSIGCONF = "CLNSIGCONF";
    public static final String BLACKLIST_BED = "BLACKLIST_BED";
    public static final String BLACKLIST_VCF = "BLACKLIST_VCF";

    public static PathogenicSummary fromContext(final VariantContext context)
    {
        String clnSig = clnSig(context);
        String clnSigConf = clnSigConf(context);
        String clinvarInfo = !clnSigConf.isEmpty() ? clnSigConf : clnSig;

        final Pathogenicity pathogenicity =
                context.getAttributeAsBoolean(BLACKLIST_BED, false) || context.getAttributeAsBoolean(BLACKLIST_VCF, false)
                        ? Pathogenicity.BENIGN_BLACKLIST
                        : Pathogenicity.fromClinvarAnnotation(clnSig, clnSigConf);

        return new PathogenicSummary(clinvarInfo, pathogenicity);
    }

    public static String clnSig(final VariantContext context)
    {
        return String.join(",", context.getAttributeAsStringList(CLNSIG, Strings.EMPTY));
    }

    public static String clnSigConf(final VariantContext context)
    {
        return String.join(",", context.getAttributeAsStringList(CLNSIGCONF, Strings.EMPTY));
    }
}
