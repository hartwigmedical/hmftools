package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.variant.VariantConsequence.consequencesToString;

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

public class ImpactVcfWriter
{
    private final VCFFileReader mHeader;
    private final VariantContextWriter mWriter;

    public ImpactVcfWriter(final String outputVCF, final String templateVCF)
    {
        mHeader = new VCFFileReader(new File(templateVCF), false);

        mWriter = new VariantContextWriterBuilder()
                .setReferenceDictionary(mHeader.getFileHeader().getSequenceDictionary())
                .setOutputFile(outputVCF)
                .build();
    }

    public final void writeHeader(final String sageVersion)
    {
        VCFHeader newHeader = new VCFHeader(mHeader.getFileHeader());
        newHeader.addMetaDataLine(new VCFHeaderLine("SageVersion", sageVersion));

        VariantTranscriptImpact.writeHeader(newHeader);
        VariantImpactSerialiser.writeHeader(newHeader);

        mWriter.writeHeader(newHeader);
    }

    public final void writeVariant(final VariantContext context, final VariantData variant, final VariantImpact variantImpact)
    {
        List<VariantTranscriptImpact> transImpacts = Lists.newArrayList();

        for(Map.Entry<String,List<VariantTransImpact>> entry : variant.getImpacts().entrySet())
        {
            final String geneName = entry.getKey();
            final List<VariantTransImpact> geneImpacts = entry.getValue();

            for(VariantTransImpact transImpact : geneImpacts)
            {
                transImpacts.add(new VariantTranscriptImpact(
                        transImpact.TransData.GeneId, geneName, transImpact.TransData.TransName,
                        consequencesToString(transImpact.consequences()), transImpact.inSpliceRegion(),
                        transImpact.hgvsCodingChange(), transImpact.hgvsProteinChange()));
            }
        }

        VariantTranscriptImpact.writeVcfData(context, transImpacts);
        VariantImpactSerialiser.writeImpactDetails(context, variantImpact);

        mWriter.add(context);
    }

    public void close()
    {
        mWriter.close();
    }
}
