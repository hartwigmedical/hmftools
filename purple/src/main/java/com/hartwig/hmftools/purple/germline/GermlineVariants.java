package com.hartwig.hmftools.purple.germline;

import static com.hartwig.hmftools.common.purple.PurpleCommon.PURPLE_GERMLINE_VCF_SUFFIX;
import static com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact.VAR_TRANS_IMPACT_ANNOATATION;
import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.variant.VariantHeader;
import com.hartwig.hmftools.purple.config.PurpleConfig;
import com.hartwig.hmftools.purple.config.ReferenceData;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;

public class GermlineVariants
{
    private final PurpleConfig mConfig;
    private final ReferenceData mReferenceData;
    private final String mVersion;

    private final List<VariantContext> mReportableVariants;

    public GermlineVariants(final PurpleConfig config, final ReferenceData referenceData, final String version)
    {
        mReferenceData = referenceData;
        mConfig = config;
        mVersion = version;

        mReportableVariants = Lists.newArrayList();
    }

    @NotNull
    public List<VariantContext> reportableVariants()
    {
        return mReportableVariants;
    }

    public void loadReportableVariants(final String germlineVcf)
    {
        if(germlineVcf.isEmpty())
            return;

        PPL_LOGGER.info("Loading germline variants from {}", germlineVcf);

        try
        {
            VCFFileReader vcfReader = new VCFFileReader(new File(germlineVcf), false);

            for(VariantContext context : vcfReader)
            {
                if(context.getAttributeAsBoolean(VariantHeader.REPORTED_FLAG, false))
                    mReportableVariants.add(context);
            }
        }
        catch(Exception e)
        {
            PPL_LOGGER.error("failed to read germline VCF from file({}): {}", germlineVcf, e.toString());
        }
    }

    public void processAndWrite(
            final String referenceId, final String tumorSample, final String germlineVcf, final PurityAdjuster purityAdjuster,
            final List<PurpleCopyNumber> copyNumbers, final Set<String> somaticReportedGenes)
    {
        mReportableVariants.clear();

        if(germlineVcf.isEmpty())
            return;

        final String outputVCF = mConfig.OutputDir + File.separator + tumorSample + PURPLE_GERMLINE_VCF_SUFFIX;

        PPL_LOGGER.info("loading germline variants from {}", germlineVcf);

        try
        {
            VCFFileReader vcfReader = new VCFFileReader(new File(germlineVcf), false);

            boolean isPaveAnnotated = vcfReader.getFileHeader().hasInfoLine(VAR_TRANS_IMPACT_ANNOATATION);

            VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(outputVCF)
                    .setOption(htsjdk.variant.variantcontext.writer.Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                    .build();

            final Consumer<VariantContext> consumer = context ->
            {
                if(context.getAttributeAsBoolean(VariantHeader.REPORTED_FLAG, false))
                {
                    mReportableVariants.add(context);
                }
                writer.add(context);
            };

            final GermlineVariantEnrichment enrichment = new GermlineVariantEnrichment(
                    mVersion, referenceId, tumorSample, mReferenceData, purityAdjuster, copyNumbers,
                    mReferenceData.GermlineHotspots, somaticReportedGenes, !isPaveAnnotated, consumer);

            writer.writeHeader(enrichment.enrichHeader(vcfReader.getFileHeader()));
            for(VariantContext context : vcfReader)
            {
                enrichment.accept(context);
            }

            enrichment.flush();
            writer.close();
        }
        catch(Exception e)
        {
            PPL_LOGGER.error(" failed to enrich germline variants: {}", e.toString());
            e.printStackTrace();
        }
    }
}
