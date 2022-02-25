package com.hartwig.hmftools.pave;

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
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VcfWriter
{
    private final VCFFileReader mHeader;
    private final VariantContextWriter mWriter;

    public static final String PON_COUNT = "PON_COUNT";
    public static final String PON_MAX = "PON_MAX";
    public static final String GNOMAD_FREQ = "GND_FREQ";

    public static final String PASS = "PASS";
    public static final String PON_FILTER = "PON";
    public static final String PON_GNOMAD_FILTER = "PONGnomad";
    public static final String PON_ARTEFACT_FILTER = "PONArtefact";

    public VcfWriter(final String outputVCF, final String templateVCF)
    {
        mHeader = new VCFFileReader(new File(templateVCF), false);

        mWriter = new VariantContextWriterBuilder()
                .setReferenceDictionary(mHeader.getFileHeader().getSequenceDictionary())
                .setOutputFile(outputVCF)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                .build();
    }

    public final void writeHeader(final String paveVersion, boolean hasGnomadFrequency, boolean hasPon)
    {
        VCFHeader newHeader = new VCFHeader(mHeader.getFileHeader());
        newHeader.addMetaDataLine(new VCFHeaderLine("PaveVersion", paveVersion));

        VariantTranscriptImpact.writeHeader(newHeader);
        VariantImpactSerialiser.writeHeader(newHeader);

        if(hasPon)
        {
            newHeader.addMetaDataLine(new VCFInfoHeaderLine(
                    PON_COUNT, 1, VCFHeaderLineType.Integer, "Cohort frequency for variant"));

            newHeader.addMetaDataLine(new VCFInfoHeaderLine(
                    PON_MAX, 1, VCFHeaderLineType.Integer, "Max read depth in any sample with variant"));

            newHeader.addMetaDataLine(new VCFFilterHeaderLine(PON_ARTEFACT_FILTER, "Filter PON artefact"));
            newHeader.addMetaDataLine(new VCFFilterHeaderLine(PON_FILTER, "Filter PON variant"));
        }

        if(hasGnomadFrequency)
        {
            newHeader.addMetaDataLine(new VCFInfoHeaderLine(
                    GNOMAD_FREQ, 1, VCFHeaderLineType.Float, "Gnomad variant frequency"));

            newHeader.addMetaDataLine(new VCFFilterHeaderLine(PON_GNOMAD_FILTER, "Filter Gnoamd PON"));
        }

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
