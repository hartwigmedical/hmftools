package com.hartwig.hmftools.common.variant.enrich;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting.ANY;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting.NONE;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting.VARIANT_NOT_LOST;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting.WILDTYPE_LOST;
import static com.hartwig.hmftools.common.variant.VariantHeader.REPORTED_DESC;
import static com.hartwig.hmftools.common.variant.VariantHeader.REPORTED_FLAG;

import java.util.Collections;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.pathogenic.Pathogenic;
import com.hartwig.hmftools.common.pathogenic.PathogenicSummary;
import com.hartwig.hmftools.common.pathogenic.PathogenicSummaryFactory;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffSummary;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class GermlineReportedEnrichment implements VariantContextEnrichment {

    private static final Set<CodingEffect> REPORTABLE_EFFECT = EnumSet.of(CodingEffect.NONSENSE_OR_FRAMESHIFT, CodingEffect.SPLICE);
    private static final double MIN_VARIANT_COPY_NUMBER = 0.5;

    @NotNull
    private final Consumer<VariantContext> consumer;
    private final Map<String, DriverGene> driverGeneMap;
    private final Set<String> somaticKnockouts;
    private final List<VariantContextDecorator> buffer = Lists.newArrayList();

    public GermlineReportedEnrichment(@NotNull final List<DriverGene> driverGenes, @NotNull final List<VariantContext> somaticVariants,
            @NotNull final Consumer<VariantContext> consumer) {
        this.consumer = consumer;

        this.driverGeneMap = driverGenes.stream().filter(DriverGene::reportGermline).collect(Collectors.toMap(DriverGene::gene, x -> x));
        this.somaticKnockouts = somaticVariants.stream()
                .map(VariantContextDecorator::new)
                .filter(VariantContextDecorator::reported)
                .map(VariantContextDecorator::gene)
                .collect(Collectors.toSet());
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
        final Map<String, Long> germlineGeneHits = buffer.stream().filter(x -> driverGeneMap.containsKey(x.gene())).filter(x -> {
            DriverGene driverGene = driverGeneMap.get(x.gene());
            return report(x,
                    downgradeWildType(driverGene.reportGermlineHotspot()),
                    downgradeWildType(driverGene.reportGermlineVariant()),
                    Collections.emptySet());
        }).map(VariantContextDecorator::gene).collect(Collectors.groupingBy(Function.identity(), Collectors.counting()));

        final Set<String> multipleGermlineGeneHits =
                germlineGeneHits.entrySet().stream().filter(x -> x.getValue() > 1).map(Map.Entry::getKey).collect(Collectors.toSet());

        final Set<String> genesWithMultipleHits = Sets.newHashSet();
        genesWithMultipleHits.addAll(multipleGermlineGeneHits);
        genesWithMultipleHits.addAll(somaticKnockouts);

        for (VariantContextDecorator variant : buffer) {
            if (report(variant, genesWithMultipleHits)) {
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

    private boolean report(@NotNull VariantContextDecorator variant, @NotNull Set<String> genesWithMultipleHits) {
        if (variant.gene().isEmpty()) {
            return false;
        }

        if (!driverGeneMap.containsKey(variant.gene())) {
            return false;
        }

        final DriverGene driverGene = driverGeneMap.get(variant.gene());
        return report(variant, driverGene.reportGermlineHotspot(), driverGene.reportGermlineVariant(), genesWithMultipleHits);
    }

    private boolean report(@NotNull VariantContextDecorator variant, @NotNull DriverGeneGermlineReporting hotspotReporting,
            @NotNull DriverGeneGermlineReporting variantReporting, @NotNull Set<String> genesWithMultipleHits) {
        if (!variant.isPass()) {
            return false;
        }

        final boolean isHotspot = variant.hotspot() == Hotspot.HOTSPOT;
        final SnpEffSummary snpEffSummary = variant.snpEffSummary();
        final PathogenicSummary pathogenicSummary = PathogenicSummaryFactory.fromContext(variant.context());
        if (!isPathogenic(isHotspot, pathogenicSummary, snpEffSummary)) {
            return false;
        }

        final DriverGeneGermlineReporting reporting = isHotspot ? hotspotReporting : variantReporting;
        if (reporting == NONE) {
            return false;
        }

        if (reporting == ANY) {
            return true;
        }

        if (variant.isVariantLost(MIN_VARIANT_COPY_NUMBER)) {
            return false;
        }

        if (reporting == VARIANT_NOT_LOST) {
            return true;
        }

        if (reporting == WILDTYPE_LOST) {
            return variant.biallelic() || genesWithMultipleHits.contains(variant.gene());
        }

        return false;
    }

    private static boolean isPathogenic(boolean isHotspot, @NotNull PathogenicSummary pathogenicSummary,
            @NotNull SnpEffSummary snpEffSummary) {
        if (pathogenicSummary.pathogenicity() == Pathogenic.BENIGN_BLACKLIST) {
            return false;
        }

        if (isHotspot || pathogenicSummary.pathogenicity().isPathogenic()) {
            return true;
        }

        return pathogenicSummary.pathogenicity() == Pathogenic.UNKNOWN && REPORTABLE_EFFECT.contains(snpEffSummary.canonicalCodingEffect());
    }
}
