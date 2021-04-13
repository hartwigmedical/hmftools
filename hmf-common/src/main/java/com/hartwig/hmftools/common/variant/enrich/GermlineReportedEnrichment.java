package com.hartwig.hmftools.common.variant.enrich;

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

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class GermlineReportedEnrichment implements VariantContextEnrichment {

    private static final double MIN_VARIANT_COPY_NUMBER = 0.5;

    @NotNull
    private final Consumer<VariantContext> consumer;
    private final Map<String, DriverGene> driverGeneMap;
    private final Set<String> somaticKnockouts;
    private final List<VariantContextDecorator> buffer = Lists.newArrayList();

    public GermlineReportedEnrichment(@NotNull final List<DriverGene> driverGenes, @NotNull final Set<String> somaticReportedGenes,
            @NotNull final Consumer<VariantContext> consumer) {
        this.consumer = consumer;

        this.driverGeneMap = driverGenes.stream().filter(DriverGene::reportGermline).collect(Collectors.toMap(DriverGene::gene, x -> x));
        this.somaticKnockouts = somaticReportedGenes;
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        buffer.add(new VariantContextDecorator(context));
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template) {
        template.addMetaDataLine(new VCFInfoHeaderLine(REPORTED_FLAG, 0, VCFHeaderLineType.Flag, REPORTED_DESC));
        return template;
    }

    public void flush() {
        final List<VariantContextDecorator> germlineHits = buffer.stream().filter(x -> driverGeneMap.containsKey(x.gene())).filter(x -> {
            DriverGene driverGene = driverGeneMap.get(x.gene());
            return report(x,
                    downgradeWildType(driverGene.reportGermlineHotspot()),
                    downgradeWildType(driverGene.reportGermlineVariant()),
                    Collections.emptySet());
        }).collect(Collectors.toList());

        for (VariantContextDecorator variant : buffer) {
            final Set<String> otherGermlineHits = germlineHits.stream()
                    .filter(x -> !x.equals(variant))
                    .filter(x -> x.gene().equals(variant.gene()))
                    .filter(x -> variant.localPhaseSet() == null || x.localPhaseSet() == null || !Objects.equals(variant.localPhaseSet(),
                            x.localPhaseSet()))
                    .map(VariantContextDecorator::gene)
                    .collect(Collectors.toSet());
            final Set<String> genesWithMultipleUnphasedHits = Sets.newHashSet();
            genesWithMultipleUnphasedHits.addAll(somaticKnockouts);
            genesWithMultipleUnphasedHits.addAll(otherGermlineHits);

            if (report(variant, genesWithMultipleUnphasedHits)) {
                variant.context().getCommonInfo().putAttribute(REPORTED_FLAG, true);
            }
            consumer.accept(variant.context());
        }

        buffer.clear();
    }

    @NotNull
    private static DriverGeneGermlineReporting downgradeWildType(@NotNull DriverGeneGermlineReporting reporting) {
        return reporting == WILDTYPE_LOST ? VARIANT_NOT_LOST : reporting;
    }

    private boolean report(@NotNull VariantContextDecorator variant, @NotNull Set<String> genesWithMultipleUnphasedHits) {
        if (variant.gene().isEmpty()) {
            return false;
        }

        if (!driverGeneMap.containsKey(variant.gene())) {
            return false;
        }

        final DriverGene driverGene = driverGeneMap.get(variant.gene());
        return report(variant, driverGene.reportGermlineHotspot(), driverGene.reportGermlineVariant(), genesWithMultipleUnphasedHits);
    }

    private boolean report(@NotNull VariantContextDecorator variant, @NotNull DriverGeneGermlineReporting hotspotReporting,
            @NotNull DriverGeneGermlineReporting variantReporting, @NotNull Set<String> genesWithMultipleUnphasedHits) {
        if (!variant.isPass()) {
            return false;
        }

        if (!variant.isPathogenic()) {
            return false;
        }

        final DriverGeneGermlineReporting reporting = variant.isHotspot() ? hotspotReporting : variantReporting;
        if (reporting == NONE) {
            return false;
        }

        if (reporting == ANY) {
            return true;
        }

        if (isVariantLost(variant, MIN_VARIANT_COPY_NUMBER)) {
            return false;
        }

        if (reporting == VARIANT_NOT_LOST) {
            return true;
        }

        if (reporting == WILDTYPE_LOST) {
            return variant.biallelic() || genesWithMultipleUnphasedHits.contains(variant.gene());
        }

        return false;
    }

    private static boolean isVariantLost(@NotNull VariantContextDecorator variant, double minVariantCopyNumber) {
        return Doubles.lessThan(variant.variantCopyNumber(), minVariantCopyNumber);
    }
}
