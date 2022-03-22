package com.hartwig.hmftools.purple.germline;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting.ANY;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting.NONE;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting.VARIANT_NOT_LOST;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting.WILDTYPE_LOST;
import static com.hartwig.hmftools.common.variant.VariantHeader.REPORTED_DESC;
import static com.hartwig.hmftools.common.variant.VariantHeader.REPORTED_FLAG;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichment;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class GermlineReportedEnrichment implements VariantContextEnrichment
{
    private static final double MIN_VARIANT_COPY_NUMBER = 0.5;

    private final Consumer<VariantContext> mConsumer;
    private final Map<String, DriverGene> mDriverGeneMap;
    private final Set<String> mSomaticKnockouts;
    private final List<VariantContextDecorator> buffer = Lists.newArrayList();

    public GermlineReportedEnrichment(
            final List<DriverGene> driverGenes, final Set<String> somaticReportedGenes, final Consumer<VariantContext> consumer)
    {
        mConsumer = consumer;
        mDriverGeneMap = driverGenes.stream().filter(DriverGene::reportGermline).collect(Collectors.toMap(DriverGene::gene, x -> x));
        mSomaticKnockouts = somaticReportedGenes;
    }

    @Override
    public void accept(final VariantContext context)
    {
        buffer.add(new VariantContextDecorator(context));
    }

    @Override
    public VCFHeader enrichHeader(final VCFHeader template)
    {
        template.addMetaDataLine(new VCFInfoHeaderLine(REPORTED_FLAG, 0, VCFHeaderLineType.Flag, REPORTED_DESC));
        return template;
    }

    public void flush()
    {
        final List<VariantContextDecorator> germlineHits = buffer.stream().filter(x -> mDriverGeneMap.containsKey(x.gene())).filter(x ->
        {
            DriverGene driverGene = mDriverGeneMap.get(x.gene());
            return report(x,
                    downgradeWildType(driverGene.reportGermlineHotspot()),
                    downgradeWildType(driverGene.reportGermlineVariant()),
                    Collections.emptySet());
        }).collect(Collectors.toList());

        for(VariantContextDecorator variant : buffer)
        {
            final Set<String> otherGermlineHits = germlineHits.stream()
                    .filter(x -> !x.equals(variant))
                    .filter(x -> x.gene().equals(variant.gene()))
                    .filter(x -> variant.localPhaseSet() == null || x.localPhaseSet() == null || !Objects.equals(variant.localPhaseSet(),
                            x.localPhaseSet()))
                    .map(VariantContextDecorator::gene)
                    .collect(Collectors.toSet());
            final Set<String> genesWithMultipleUnphasedHits = Sets.newHashSet();
            genesWithMultipleUnphasedHits.addAll(mSomaticKnockouts);
            genesWithMultipleUnphasedHits.addAll(otherGermlineHits);

            if(report(variant, genesWithMultipleUnphasedHits))
            {
                variant.context().getCommonInfo().putAttribute(REPORTED_FLAG, true);
            }
            mConsumer.accept(variant.context());
        }

        buffer.clear();
    }

    @NotNull
    private static DriverGeneGermlineReporting downgradeWildType(DriverGeneGermlineReporting reporting)
    {
        return reporting == WILDTYPE_LOST ? VARIANT_NOT_LOST : reporting;
    }

    private boolean report(VariantContextDecorator variant, Set<String> genesWithMultipleUnphasedHits)
    {
        if(variant.gene().isEmpty())
        {
            return false;
        }

        if(!mDriverGeneMap.containsKey(variant.gene()))
        {
            return false;
        }

        final DriverGene driverGene = mDriverGeneMap.get(variant.gene());
        return report(variant, driverGene.reportGermlineHotspot(), driverGene.reportGermlineVariant(), genesWithMultipleUnphasedHits);
    }

    private boolean report(VariantContextDecorator variant, DriverGeneGermlineReporting hotspotReporting,
            DriverGeneGermlineReporting variantReporting, Set<String> genesWithMultipleUnphasedHits)
    {
        if(!variant.isPass())
        {
            return false;
        }

        if(!variant.isPathogenic())
        {
            return false;
        }

        final DriverGeneGermlineReporting reporting = variant.isHotspot() ? hotspotReporting : variantReporting;
        if(reporting == NONE)
        {
            return false;
        }

        if(reporting == ANY)
        {
            return true;
        }

        if(isVariantLost(variant, MIN_VARIANT_COPY_NUMBER))
        {
            return false;
        }

        if(reporting == VARIANT_NOT_LOST)
        {
            return true;
        }

        if(reporting == WILDTYPE_LOST)
        {
            return variant.biallelic() || genesWithMultipleUnphasedHits.contains(variant.gene());
        }

        return false;
    }

    private static boolean isVariantLost(VariantContextDecorator variant, double minVariantCopyNumber)
    {
        return Doubles.lessThan(variant.variantCopyNumber(), minVariantCopyNumber);
    }
}
