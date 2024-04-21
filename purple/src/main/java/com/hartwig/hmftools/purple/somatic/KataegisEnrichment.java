package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.KATAEGIS_FLAG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.KATAEGIS_FLAG_DESCRIPTION;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TRINUCLEOTIDE_CONTEXT;

import java.util.concurrent.atomic.AtomicInteger;

import org.apache.logging.log4j.util.Strings;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class KataegisEnrichment
{
    private final KataegisQueue mForwardDetector;
    private final KataegisQueue mReverseDetector;

    public KataegisEnrichment(final AtomicInteger kataegisId)
    {
        mReverseDetector = new KataegisQueue("REV", kataegisId, KataegisEnrichment::isReverseCandidate, null);
        mForwardDetector = new KataegisQueue("FWD", kataegisId, KataegisEnrichment::isForwardCandidate, mReverseDetector::processVariant);
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

    public static void enrichHeader(final VCFHeader template)
    {
        template.addMetaDataLine(new VCFInfoHeaderLine(KATAEGIS_FLAG, 1, VCFHeaderLineType.String, KATAEGIS_FLAG_DESCRIPTION));
    }

    private static boolean isForwardCandidate(final SomaticVariant variant)
    {
        final VariantContext context = variant.context();

        final boolean altMatch = context.getAlternateAlleles().stream()
                .anyMatch(x -> x.getBaseString().equals("T") || x.getBaseString().equals("G"));

        final String triContext = context.getAttributeAsString(TRINUCLEOTIDE_CONTEXT, Strings.EMPTY);
        final boolean triMatch = triContext.startsWith("TC");

        return variant.isPass() && triMatch && altMatch;
    }

    private static boolean isReverseCandidate(final SomaticVariant variant)
    {
        final VariantContext context = variant.context();

        final boolean altMatch = context.getAlternateAlleles().stream()
                .anyMatch(x -> x.getBaseString().equals("C") || x.getBaseString().equals("A"));

        final String triContext = context.getAttributeAsString(TRINUCLEOTIDE_CONTEXT, Strings.EMPTY);
        final boolean triMatch = triContext.endsWith("GA");

        return variant.isPass() && triMatch && altMatch;
    }
}
