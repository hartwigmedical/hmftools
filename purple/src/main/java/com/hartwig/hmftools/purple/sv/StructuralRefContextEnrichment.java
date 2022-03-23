package com.hartwig.hmftools.purple.sv;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.REF_CONTEXT_FLAG;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.relativePositionAndRef;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichment;

import org.apache.commons.math3.util.Pair;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class StructuralRefContextEnrichment implements VariantContextEnrichment
{
    private static final int REF_CONTEXT_DISTANCE = 10;

    private static final String REF_CONTEXT_DESCRIPTION = "Reference genome surrounding break";

    private final IndexedFastaSequenceFile reference;
    private final Consumer<VariantContext> consumer;

    public StructuralRefContextEnrichment(final IndexedFastaSequenceFile reference, final Consumer<VariantContext> consumer)
    {
        this.reference = reference;
        this.consumer = consumer;
    }

    @Override
    public VCFHeader enrichHeader(final VCFHeader template)
    {
        template.addMetaDataLine(new VCFInfoHeaderLine(REF_CONTEXT_FLAG, 1, VCFHeaderLineType.String, REF_CONTEXT_DESCRIPTION));

        return template;
    }

    @Override
    public void accept(final VariantContext context)
    {
        final Pair<Integer, String> relativePositionAndRef = relativePositionAndRef(reference, context);
        if(relativePositionAndRef.getFirst() > -1)
        {
            addRefContext(context, relativePositionAndRef);
        }
        consumer.accept(context);
    }

    @Override
    public void flush() { }

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
