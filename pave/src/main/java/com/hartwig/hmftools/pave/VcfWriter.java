package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.pave.external.GnomadAnnotation.GNOMAD_VCF_TAG;

import java.io.File;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;
import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VcfWriter
{
    private final VCFFileReader mHeader;
    private final VariantContextWriter mWriter;

    public VcfWriter(final String outputVCF, final String templateVCF)
    {
        mHeader = new VCFFileReader(new File(templateVCF), false);

        mWriter = new VariantContextWriterBuilder()
                .setReferenceDictionary(mHeader.getFileHeader().getSequenceDictionary())
                .setOutputFile(outputVCF)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                .build();
    }

    public final void writeHeader(final String paveVersion, boolean hasGnomadFrequency)
    {
        VCFHeader newHeader = new VCFHeader(mHeader.getFileHeader());
        newHeader.addMetaDataLine(new VCFHeaderLine("PaveVersion", paveVersion));

        VariantTranscriptImpact.writeHeader(newHeader);
        VariantImpactSerialiser.writeHeader(newHeader);

        if(hasGnomadFrequency)
        {
            newHeader.addMetaDataLine(new VCFInfoHeaderLine(
                    GNOMAD_VCF_TAG, 1, VCFHeaderLineType.Float, "Gnomad variant frequency"));
        }

        mWriter.writeHeader(newHeader);
    }

    public final void writeVariant(final VariantContext context, final VariantData variant, final VariantImpact variantImpact)
    {
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

            VariantTranscriptImpact.writeVcfData(context, transImpacts);

            if(variantImpact != null)
                VariantImpactSerialiser.writeImpactDetails(context, variantImpact);

            if(variant.gnomadFrequency() != null)
                context.getCommonInfo().putAttribute(GNOMAD_VCF_TAG, variant.gnomadFrequency());
        }

        mWriter.add(context);
    }

    public void close()
    {
        mWriter.close();
    }
}
