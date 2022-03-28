package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.KATAEGIS_FLAG;
import static com.hartwig.hmftools.common.variant.VariantHeader.PASS;

import java.util.Set;

import com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment;

import org.apache.logging.log4j.util.Strings;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class KataegisEnrichment
{
    private static final String KATAEGIS_FLAG_DESCRIPTION = "Forward/reverse kataegis id";

    private final KataegisQueue mForwardDetector;
    private final KataegisQueue mReverseDetector;

    public KataegisEnrichment()
    {
        mReverseDetector = new KataegisQueue("REV", KataegisEnrichment::isReverseCandidate, null);
        mForwardDetector = new KataegisQueue("FWD", KataegisEnrichment::isForwardCandidate, mReverseDetector::processVariant);
    }

    public void processVariant(final SomaticVariant variant)
    {
        mForwardDetector.processVariant(variant);
    }

    public void flush()
    {
        mForwardDetector.flush();
        mReverseDetector.flush();
    }

    public static VCFHeader enrichHeader(final VCFHeader template)
    {
        template.addMetaDataLine(new VCFInfoHeaderLine(KATAEGIS_FLAG, 1, VCFHeaderLineType.String, KATAEGIS_FLAG_DESCRIPTION));
        return template;
    }

    private static boolean isForwardCandidate(final SomaticVariant variant)
    {
        final VariantContext context = variant.context();

        final boolean altMatch = context.getAlternateAlleles().stream()
                .anyMatch(x -> x.getBaseString().equals("T") || x.getBaseString().equals("G"));

        final String triContext = context.getAttributeAsString(SomaticRefContextEnrichment.TRINUCLEOTIDE_FLAG, Strings.EMPTY);
        final boolean triMatch = triContext.startsWith("TC");

        return variant.isPass() && triMatch && altMatch;
    }

    private static boolean isReverseCandidate(final SomaticVariant variant)
    {
        final VariantContext context = variant.context();

        final boolean altMatch = context.getAlternateAlleles().stream()
                .anyMatch(x -> x.getBaseString().equals("C") || x.getBaseString().equals("A"));

        final String triContext = context.getAttributeAsString(SomaticRefContextEnrichment.TRINUCLEOTIDE_FLAG, Strings.EMPTY);
        final boolean triMatch = triContext.endsWith("GA");

        return variant.isPass() && triMatch && altMatch;
    }
}
