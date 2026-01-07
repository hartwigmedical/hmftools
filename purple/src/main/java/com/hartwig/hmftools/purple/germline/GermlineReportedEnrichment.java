package com.hartwig.hmftools.purple.germline;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_DESC;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_FLAG;
import static com.hartwig.hmftools.purple.drivers.SomaticVariantDrivers.addReportableTranscriptList;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class GermlineReportedEnrichment
{
    private final Map<String,DriverGene> mDriverGeneMap;
    private final List<GermlineVariant> mBuffer;

    public GermlineReportedEnrichment(final Map<String,DriverGene> driverGenes)
    {
        mDriverGeneMap = driverGenes;
        mBuffer = Lists.newArrayList();
    }

    public void processVariant(final GermlineVariant variant)
    {
        mBuffer.add(variant);
    }

    public static void enrichHeader(final VCFHeader header)
    {
        header.addMetaDataLine(new VCFInfoHeaderLine(REPORTED_FLAG, 0, VCFHeaderLineType.Flag, REPORTED_DESC));
    }

    public void flush()
    {
        for(GermlineVariant variant : mBuffer)
        {
            if(!mDriverGeneMap.containsKey(variant.gene()))
                continue;

            if(isCandidateReportable(variant.decorator()))
            {
                variant.context().getCommonInfo().putAttribute(REPORTED_FLAG, true);
                addReportableTranscriptList(variant.type(), variant.context(), variant.variantImpact());
            }
        }

        mBuffer.clear();
    }

    public boolean isCandidateReportable(final VariantContextDecorator variant)
    {
        if(variant.gene().isEmpty())
            return false;

        if(!mDriverGeneMap.containsKey(variant.gene()))
            return false;

        if(!variant.isPass())
            return false;

        return variant.isGermlinePathogenic();
    }

    /*
    public static boolean isReportable(
            final GermlineVariant variant, final DriverGeneGermlineReporting hotspotReporting,
            final DriverGeneGermlineReporting variantReporting, final Set<String> genesWithMultipleUnphasedHits)
    {
        DriverGeneGermlineReporting reporting = variant.isHotspot() ? hotspotReporting : variantReporting;

        if(reporting == NONE)
            return false;

        if(reporting == ANY)
            return true;

        if(isVariantLost(variant, MIN_VARIANT_COPY_NUMBER))
            return false;

        if(reporting == VARIANT_NOT_LOST)
            return true;

        if(reporting == WILDTYPE_LOST)
            return variant.biallelic() || genesWithMultipleUnphasedHits.contains(variant.gene());

        return false;
    }
    */

    private static boolean isVariantLost(final GermlineVariant variant, double minVariantCopyNumber)
    {
        return Doubles.lessThan(variant.decorator().variantCopyNumber(), minVariantCopyNumber);
    }
}
