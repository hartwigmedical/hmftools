package com.hartwig.hmftools.purple.germline;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting.ANY;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting.NONE;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting.VARIANT_NOT_LOST;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting.WILDTYPE_LOST;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_DESC;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_FLAG;
import static com.hartwig.hmftools.purple.drivers.SomaticVariantDrivers.addReportableTranscriptList;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class GermlineReportedEnrichment
{
    private static final double MIN_VARIANT_COPY_NUMBER = 0.5;

    private final Map<String, DriverGene> mDriverGeneMap;
    private final Set<String> mSomaticKnockouts;
    private final List<GermlineVariant> mBuffer;

    public GermlineReportedEnrichment(final List<DriverGene> driverGenes, final Set<String> somaticReportedGenes)
    {
        mDriverGeneMap = driverGenes.stream().filter(DriverGene::reportGermline).collect(Collectors.toMap(DriverGene::gene, x -> x));
        mSomaticKnockouts = somaticReportedGenes;
        mBuffer = Lists.newArrayList();
    }

    public void processVariant(final GermlineVariant variant)
    {
        // critical that the new context and filters are used
        mBuffer.add(variant);
    }

    public static void enrichHeader(final VCFHeader header)
    {
        header.addMetaDataLine(new VCFInfoHeaderLine(REPORTED_FLAG, 0, VCFHeaderLineType.Flag, REPORTED_DESC));
    }

    public void flush()
    {
        final List<GermlineVariant> germlineHits = Lists.newArrayList();

        for(GermlineVariant variant : mBuffer)
        {
            if(!mDriverGeneMap.containsKey(variant.gene()))
                continue;

            DriverGene driverGene = mDriverGeneMap.get(variant.gene());

            if(isReportable(
                    variant.decorator(), downgradeWildType(driverGene.reportGermlineHotspot()),
                    downgradeWildType(driverGene.reportGermlineVariant()), Collections.emptySet()))
            {
                germlineHits.add(variant);
            }
        }

        for(GermlineVariant variant : mBuffer)
        {
            final Set<String> otherGermlineHits = germlineHits.stream()
                    .filter(x -> !x.equals(variant))
                    .filter(x -> x.gene().equals(variant.gene()))
                    .filter(x -> variant.decorator().localPhaseSet() == null || x.decorator().localPhaseSet() == null
                            || !Objects.equals(variant.decorator().localPhaseSet(), x.decorator().localPhaseSet()))
                    .map(x -> x.gene())
                    .collect(Collectors.toSet());

            final Set<String> genesWithMultipleUnphasedHits = Sets.newHashSet();
            genesWithMultipleUnphasedHits.addAll(mSomaticKnockouts);
            genesWithMultipleUnphasedHits.addAll(otherGermlineHits);

            if(report(variant.decorator(), genesWithMultipleUnphasedHits))
            {
                variant.context().getCommonInfo().putAttribute(REPORTED_FLAG, true);

                addReportableTranscriptList(variant.type(), variant.context(), variant.variantImpact());
            }
        }

        mBuffer.clear();
    }

    private static DriverGeneGermlineReporting downgradeWildType(DriverGeneGermlineReporting reporting)
    {
        return reporting == WILDTYPE_LOST ? VARIANT_NOT_LOST : reporting;
    }

    public boolean report(final VariantContextDecorator variant, final Set<String> genesWithMultipleUnphasedHits)
    {
        if(variant.gene().isEmpty())
            return false;

        if(!mDriverGeneMap.containsKey(variant.gene()))
            return false;

        final DriverGene driverGene = mDriverGeneMap.get(variant.gene());

        return isReportable(variant, driverGene.reportGermlineHotspot(), driverGene.reportGermlineVariant(), genesWithMultipleUnphasedHits);
    }

    public static boolean isReportable(
            final VariantContextDecorator variant, final DriverGeneGermlineReporting hotspotReporting,
            final DriverGeneGermlineReporting variantReporting, final Set<String> genesWithMultipleUnphasedHits)
    {
        if(!variant.isPass())
            return false;

        if(!variant.isPathogenic())
            return false;

        final DriverGeneGermlineReporting reporting = variant.isHotspot() ? hotspotReporting : variantReporting;

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

    private static boolean isVariantLost(VariantContextDecorator variant, double minVariantCopyNumber)
    {
        return Doubles.lessThan(variant.variantCopyNumber(), minVariantCopyNumber);
    }
}
