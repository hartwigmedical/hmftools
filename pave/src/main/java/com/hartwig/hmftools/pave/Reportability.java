package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_DESC;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_FLAG;
import static com.hartwig.hmftools.common.variant.Hotspot.HOTSPOT_DESCRIPTION;
import static com.hartwig.hmftools.common.variant.Hotspot.HOTSPOT_FLAG;

import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.ReportablePredicate;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class Reportability
{
    private final ReportablePredicate mReportableOncoGenes;
    private final ReportablePredicate mReportableTsgGenes;

    public Reportability(final List<DriverGene> driverGenes)
    {
        mReportableOncoGenes = new ReportablePredicate(DriverCategory.ONCO, driverGenes);
        mReportableTsgGenes = new ReportablePredicate(DriverCategory.TSG, driverGenes);
    }

    public boolean isReported(final VariantData variant, final VariantImpact variantImpact, boolean isHotspot)
    {
        if(variantImpact == null)
            return false;

        if(mReportableOncoGenes.isReportable(variantImpact, variant.type(), variant.repeatCount(), isHotspot))
        {
            return true;
        }

        if(mReportableTsgGenes.isReportable(variantImpact, variant.type(), variant.repeatCount(), isHotspot))
        {
            return true;
        }

        return false;
    }

    public void setReportability(final VariantData variant, final VariantImpact variantImpact)
    {
        boolean isHotspot = variant.tier() == VariantTier.HOTSPOT;

        if(isHotspot)
        {
            variant.context().getCommonInfo().putAttribute(HOTSPOT_FLAG, true);
        }

        if(isReported(variant, variantImpact, isHotspot))
        {
            variant.context().getCommonInfo().putAttribute(REPORTED_FLAG, true);
        }
    }

    public static void addHeader(final VCFHeader header)
    {
        header.addMetaDataLine(new VCFInfoHeaderLine(HOTSPOT_FLAG, 0, VCFHeaderLineType.Flag, HOTSPOT_DESCRIPTION));
        header.addMetaDataLine(new VCFInfoHeaderLine(REPORTED_FLAG, 0, VCFHeaderLineType.Flag, REPORTED_DESC));
   }
}
