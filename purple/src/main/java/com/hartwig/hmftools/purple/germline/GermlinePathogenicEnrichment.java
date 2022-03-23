package com.hartwig.hmftools.purple.germline;

import java.util.Arrays;
import java.util.StringJoiner;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.pathogenic.PathogenicSummary;
import com.hartwig.hmftools.common.pathogenic.PathogenicSummaryFactory;
import com.hartwig.hmftools.common.pathogenic.Pathogenicity;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class GermlinePathogenicEnrichment
{
    private static final String PATHOGENICITY = "PATH";

    public static void processVariant(final VariantContext context)
    {
        final PathogenicSummary summary = PathogenicSummaryFactory.fromContext(context);
        context.getCommonInfo().putAttribute(PATHOGENICITY, summary.pathogenicity().toString());
    }

    public static VCFHeader enrichHeader(final VCFHeader template)
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
