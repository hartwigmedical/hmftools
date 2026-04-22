package com.hartwig.hmftools.common.pathogenic;

import static com.hartwig.hmftools.common.variant.PaveVcfTags.BLACKLIST_BED_FLAG;
import static com.hartwig.hmftools.common.variant.PaveVcfTags.BLACKLIST_VCF_FLAG;

import org.apache.logging.log4j.util.Strings;

import htsjdk.variant.variantcontext.VariantContext;

public final class PathogenicSummaryFactory
{
    public static final String CLNSIG = "CLNSIG";
    public static final String CLNSIGCONF = "CLNSIGCONF";

    public static PathogenicSummary fromContext(final VariantContext context)
    {
        String clnSig = clnSig(context);
        String clnSigConf = clnSigConf(context);
        String clinvarInfo = !clnSigConf.isEmpty() ? clnSigConf : clnSig;

        boolean isBlacklisted = context.getAttributeAsBoolean(BLACKLIST_BED_FLAG, false)
                || context.getAttributeAsBoolean(BLACKLIST_VCF_FLAG, false);

        final Pathogenicity pathogenicity = isBlacklisted ?
                Pathogenicity.BENIGN_BLACKLIST : Pathogenicity.fromClinvarAnnotation(clnSig, clnSigConf);

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
