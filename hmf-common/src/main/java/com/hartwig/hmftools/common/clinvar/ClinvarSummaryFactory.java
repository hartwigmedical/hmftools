package com.hartwig.hmftools.common.clinvar;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class ClinvarSummaryFactory {

    public static final String CLNSIG = "CLNSIG";
    public static final String CLNSIGCONF = "CLNSIGCONF";

    @NotNull
    public static ClinvarSummary fromContext(@NotNull VariantContext context) {

        final String clnSig = context.getAttributeAsString(CLNSIG, Strings.EMPTY);
        final String clnSigConf = String.join(",", context.getAttributeAsStringList(CLNSIGCONF, Strings.EMPTY));
        final ClinvarPathogenicity clnPath = ClinvarPathogenicity.fromClinvarAnnotation(clnSig, clnSigConf);
        final String clinvarInfo = !clnSigConf.isEmpty() ? clnSigConf : clnSig;

        return ImmutableClinvarSummary.builder().info(clinvarInfo).pathogenicity(clnPath).build();
    }

}
