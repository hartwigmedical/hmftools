package com.hartwig.hmftools.common.variant.enrich;

import java.util.Arrays;
import java.util.StringJoiner;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.clinvar.ClinvarPathogenicity;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class ClinvarEnrichment implements VariantContextEnrichment  {

    private static final String CLNSIG = "CLNSIG";
    private static final String CLNSIGCONF = "CLNSIGCONF";
    private static final String CLINVAR_PATHOGENICITY = "CLNPATH";

    private final Consumer<VariantContext> consumer;

    public ClinvarEnrichment(final Consumer<VariantContext> consumer) {
        this.consumer = consumer;
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        final String clnSig = context.getAttributeAsString(CLNSIG, Strings.EMPTY);
        final String clnSigConf = context.getAttributeAsString(CLNSIGCONF, Strings.EMPTY);
        final String clnPath = ClinvarPathogenicity.fromClinvarAnnotation(clnSig, clnSigConf).toString();

        context.getCommonInfo().putAttribute(CLINVAR_PATHOGENICITY, clnPath);
        consumer.accept(context);
    }

    @Override
    public void flush() {

    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template) {
        StringJoiner joiner = new StringJoiner(",");
        Arrays.stream(ClinvarPathogenicity.values()).forEach(x -> joiner.add(x.toString()));
        template.addMetaDataLine(new VCFInfoHeaderLine(CLINVAR_PATHOGENICITY, 1, VCFHeaderLineType.String, "Clinical pathogenicity [" + joiner.toString() + "]"));
        return template;
    }
}
