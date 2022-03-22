package com.hartwig.hmftools.purple.germline;

import java.util.Arrays;
import java.util.StringJoiner;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.pathogenic.PathogenicSummary;
import com.hartwig.hmftools.common.pathogenic.PathogenicSummaryFactory;
import com.hartwig.hmftools.common.pathogenic.Pathogenicity;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichment;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class GermlinePathogenicEnrichment implements VariantContextEnrichment
{
    private static final String PATHOGENICITY = "PATH";

    private final Consumer<VariantContext> mConsumer;

    public GermlinePathogenicEnrichment(final Consumer<VariantContext> consumer)
    {
        this.mConsumer = consumer;
    }

    @Override
    public void accept(final VariantContext context)
    {
        final PathogenicSummary summary = PathogenicSummaryFactory.fromContext(context);
        context.getCommonInfo().putAttribute(PATHOGENICITY, summary.pathogenicity().toString());
        mConsumer.accept(context);
    }

    @Override
    public void flush() { }

    @Override
    public VCFHeader enrichHeader(final VCFHeader template)
    {
        StringJoiner joiner = new StringJoiner(",");
        Arrays.stream(Pathogenicity.values()).forEach(x -> joiner.add(x.toString()));
        template.addMetaDataLine(new VCFInfoHeaderLine(PATHOGENICITY,
                1,
                VCFHeaderLineType.String,
                "Pathogenicity [" + joiner.toString() + "]"));
        return template;
    }
}
