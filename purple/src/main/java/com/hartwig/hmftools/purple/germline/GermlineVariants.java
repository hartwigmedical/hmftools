package com.hartwig.hmftools.purple.germline;

import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.File;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.purple.fitting.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.variant.CommonVcfTags;
import com.hartwig.hmftools.purple.PurpleConfig;
import com.hartwig.hmftools.purple.ReferenceData;

import com.google.common.collect.Lists;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public class GermlineVariants
{
    private final PurpleConfig mConfig;
    private final ReferenceData mReferenceData;
    private final String mVersion;

    private final List<GermlineVariant> mVariants;
    private final List<GermlineVariant> mReportableVariants;

    public GermlineVariants(final PurpleConfig config, final ReferenceData referenceData, final String version)
    {
        mReferenceData = referenceData;
        mConfig = config;
        mVersion = version;

        mVariants = Lists.newArrayList();
        mReportableVariants = Lists.newArrayList();
    }

    public List<GermlineVariant> reportableVariants()
    {
        return mReportableVariants;
    }

    private void loadGermlineVariants(final String germlineVcf)
    {
        if(germlineVcf.isEmpty())
            return;

        VcfFileReader vcfReader = new VcfFileReader(germlineVcf);

        for(VariantContext context : vcfReader.iterator())
        {
            GermlineVariant variant = new GermlineVariant(context);

            mVariants.add(variant);
        }

        PPL_LOGGER.info("load {} germline variants from {}", mVariants.size(), germlineVcf);
    }

    public void processAndWrite(
            final String referenceId, final String tumorSample, final String germlineVcf, @Nullable final PurityAdjuster purityAdjuster,
            final List<PurpleCopyNumber> copyNumbers, final Set<String> somaticReportedGenes)
    {
        mReportableVariants.clear();

        if(germlineVcf.isEmpty())
            return;

        String outputVCF = PurpleCommon.purpleGermlineVcfFile(mConfig.OutputDir, tumorSample);

        loadGermlineVariants(germlineVcf);

        VCFFileReader vcfReader = new VCFFileReader(new File(germlineVcf), false);

        VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(outputVCF)
                .setOption(htsjdk.variant.variantcontext.writer.Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                .build();

        GermlineVariantEnrichment enrichment = new GermlineVariantEnrichment(
                mVersion, referenceId, tumorSample, mReferenceData, purityAdjuster, copyNumbers,
                mReferenceData.GermlineHotspots, somaticReportedGenes);

        VCFHeader header = vcfReader.getFileHeader();
        enrichment.enrichHeader(header);
        writer.writeHeader(header);

        for(GermlineVariant variant : mVariants)
        {
            enrichment.enrichVariant(variant);
        }

        enrichment.flush();

        for(GermlineVariant variant : mVariants)
        {
            VariantContext newContext = new VariantContextBuilder(variant.context()).filters(variant.filters()).make();

            if(newContext.getAttributeAsBoolean(CommonVcfTags.REPORTED_FLAG, false))
                mReportableVariants.add(variant);

            writer.add(newContext);
        }

        writer.close();
    }
}
