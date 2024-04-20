package com.hartwig.hmftools.purple.sv;

import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_CONTEXT_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_CONTEXT_FLAG;
import static com.hartwig.hmftools.common.variant.VariantUtils.relativePositionAndRef;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.util.function.Consumer;

import org.apache.commons.math3.util.Pair;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class StructuralRefContextEnrichment
{
    private static final int REF_CONTEXT_DISTANCE = 10;

    private final IndexedFastaSequenceFile mRefGenome;
    private final Consumer<VariantContext> mConsumer;

    public StructuralRefContextEnrichment(final IndexedFastaSequenceFile reference, final Consumer<VariantContext> consumer)
    {
        mRefGenome = reference;
        mConsumer = consumer;
    }

    public void enrichHeader(final VCFHeader template)
    {
        template.addMetaDataLine(new VCFInfoHeaderLine(REF_CONTEXT_FLAG, 1, VCFHeaderLineType.String, REF_CONTEXT_DESC));
    }

    public void accept(final VariantContext context)
    {
        if(context.hasAttribute(REF_CONTEXT_FLAG))
        {
            PPL_LOGGER.debug("variant({}:{}) already has ref-context set", context.getContig(), context.getStart());
            mConsumer.accept(context);
            return;
        }

        final Pair<Integer, String> relativePositionAndRef = relativePositionAndRef(mRefGenome, context);
        if(relativePositionAndRef.getFirst() > -1)
        {
            addRefContext(context, relativePositionAndRef);
        }
        mConsumer.accept(context);
    }

    private void addRefContext(final VariantContext variant, final Pair<Integer, String> relativePositionAndRef)
    {
        final int relativePosition = relativePositionAndRef.getFirst();
        final String sequence = relativePositionAndRef.getSecond();
        if(!sequence.isEmpty())
        {
            final String tri = sequence.substring(Math.max(0, relativePosition - REF_CONTEXT_DISTANCE),
                    Math.min(sequence.length(), relativePosition + REF_CONTEXT_DISTANCE + 1));
            variant.getCommonInfo().putAttribute(REF_CONTEXT_FLAG, tri);
        }
    }
}
