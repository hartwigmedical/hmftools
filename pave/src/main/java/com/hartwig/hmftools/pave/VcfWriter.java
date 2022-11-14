package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.variant.PaveVcfTags.GNOMAD_FREQ;
import static com.hartwig.hmftools.pave.PonAnnotation.PON_COUNT;
import static com.hartwig.hmftools.pave.PonAnnotation.PON_MAX;

import java.io.File;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;
import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

public class VcfWriter
{
    private final VCFFileReader mHeader;
    private final VariantContextWriter mWriter;

    public static final String PASS = "PASS";

    public VcfWriter(final String outputVCF, final String templateVCF)
    {
        mHeader = new VCFFileReader(new File(templateVCF), false);

        mWriter = new VariantContextWriterBuilder()
                .setReferenceDictionary(mHeader.getFileHeader().getSequenceDictionary())
                .setOutputFile(outputVCF)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                .build();
    }

    public final void writeHeader(
            final String paveVersion, boolean hasGnomadFrequency, boolean hasPon, boolean hasMappability, boolean hasClinvar,
            boolean hasBlacklistings, boolean hasReportability)
    {
        VCFHeader newHeader = new VCFHeader(mHeader.getFileHeader());
        newHeader.addMetaDataLine(new VCFHeaderLine("PaveVersion", paveVersion));

        VariantTranscriptImpact.writeHeader(newHeader);
        VariantImpactSerialiser.writeHeader(newHeader);

        if(hasPon)
        {
            PonAnnotation.addHeader(newHeader);
        }

        if(hasGnomadFrequency)
        {
            GnomadAnnotation.addHeader(newHeader);
        }

        if(hasMappability)
        {
            Mappability.addHeader(newHeader);
        }

        if(hasClinvar)
        {
            ClinvarAnnotation.addHeader(newHeader);
        }

        if(hasBlacklistings)
        {
            Blacklistings.addHeader(newHeader);
        }

        if(hasReportability)
            Reportability.addHeader(newHeader);

        mWriter.writeHeader(newHeader);
    }

    public final void writeVariant(final VariantContext context, final VariantData variant, final VariantImpact variantImpact)
    {
        VariantContextBuilder builder = new VariantContextBuilder(variant.context())
                .genotypes(variant.context().getGenotypes())
                .filters(context.getFilters());

        if(!variant.filters().isEmpty())
        {
            builder.getFilters().addAll(variant.filters());
            builder.getFilters().remove(PASS);
        }

        if(builder.getFilters().isEmpty())
            builder.getFilters().add(PASS);

        final VariantContext newContext = builder.make();

        if(!variant.getImpacts().isEmpty())
        {
            List<VariantTranscriptImpact> transImpacts = Lists.newArrayList();

            for(Map.Entry<String, List<VariantTransImpact>> entry : variant.getImpacts().entrySet())
            {
                final String geneName = entry.getKey();
                final List<VariantTransImpact> geneImpacts = entry.getValue();

                for(VariantTransImpact transImpact : geneImpacts)
                {
                    transImpacts.add(new VariantTranscriptImpact(
                            transImpact.TransData.GeneId, geneName, transImpact.TransData.TransName,
                            transImpact.effectsStr(), transImpact.inSpliceRegion(),
                            transImpact.hgvsCoding(), transImpact.hgvsProtein()));
                }
            }

            VariantTranscriptImpact.writeVcfData(newContext, transImpacts);

            if(variantImpact != null)
                VariantImpactSerialiser.writeImpactDetails(newContext, variantImpact);
        }

        if(variant.gnomadFrequency() != null && !newContext.getCommonInfo().hasAttribute(GNOMAD_FREQ))
            newContext.getCommonInfo().putAttribute(GNOMAD_FREQ, variant.gnomadFrequency());

        if((variant.ponSampleCount() > 0 || variant.ponMaxReadCount() > 0) && !newContext.getCommonInfo().hasAttribute(PON_COUNT))
        {
            newContext.getCommonInfo().putAttribute(PON_COUNT, variant.ponSampleCount());
            newContext.getCommonInfo().putAttribute(PON_MAX, variant.ponMaxReadCount());
        }

        mWriter.add(newContext);
    }

    public void close()
    {
        mWriter.close();
    }
}
