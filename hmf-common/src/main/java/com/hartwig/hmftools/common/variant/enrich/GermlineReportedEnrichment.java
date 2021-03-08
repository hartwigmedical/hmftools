package com.hartwig.hmftools.common.variant.enrich;

import static com.hartwig.hmftools.common.variant.VariantHeader.REPORTED_DESC;
import static com.hartwig.hmftools.common.variant.VariantHeader.REPORTED_FLAG;

import java.util.EnumSet;
import java.util.List;
import java.util.Set;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.pathogenic.Pathogenic;
import com.hartwig.hmftools.common.pathogenic.PathogenicSummary;
import com.hartwig.hmftools.common.pathogenic.PathogenicSummaryFactory;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffSummary;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffSummaryFactory;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class GermlineReportedEnrichment implements VariantContextEnrichment {

    private static final Set<CodingEffect> REPORTABLE_EFFECT = EnumSet.of(CodingEffect.NONSENSE_OR_FRAMESHIFT, CodingEffect.SPLICE);

    @NotNull
    private final Consumer<VariantContext> consumer;
    @NotNull
    private final Set<String> reportableVariantGenes;
    @NotNull
    private final Set<String> reportableHotspotGenes;

    public GermlineReportedEnrichment(@NotNull final List<DriverGene> driverGenes, @NotNull final Consumer<VariantContext> consumer) {
        this.consumer = consumer;
        this.reportableVariantGenes =
                driverGenes.stream().filter(DriverGene::reportGermlineVariant).map(DriverGene::gene).collect(Collectors.toSet());
        this.reportableHotspotGenes =
                driverGenes.stream().filter(DriverGene::reportGermlineHotspot).map(DriverGene::gene).collect(Collectors.toSet());
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        if (report(context)) {
            context.getCommonInfo().putAttribute(REPORTED_FLAG, true);
        }
        consumer.accept(context);
    }

    private boolean report(@NotNull final VariantContext context) {
        final boolean isPass = !context.isFiltered() || (context.getFilters().size() == 1 && context.getFilters().contains("PASS"));

        if (!isPass) {
            return false;
        }

        final PathogenicSummary pathogenicSummary = PathogenicSummaryFactory.fromContext(context);
        final SnpEffSummary snpEffSummary = SnpEffSummaryFactory.fromSnpEffEnrichment(context);
        final String gene = snpEffSummary.gene();
        final boolean inHotspotGenes = reportableHotspotGenes.contains(gene);
        final boolean isHotspot = HotspotEnrichment.fromVariant(context).equals(Hotspot.HOTSPOT);
        if (isHotspot && inHotspotGenes && !pathogenicSummary.pathogenicity().equals(Pathogenic.BENIGN_BLACKLIST)) {
            return true;
        }

        final boolean inVariantGenes = reportableVariantGenes.contains(gene);
        if (inVariantGenes) {
            if (pathogenicSummary.pathogenicity().isPathogenic()) {
                return true;
            }

            if (pathogenicSummary.pathogenicity().equals(Pathogenic.UNKNOWN)) {
                return REPORTABLE_EFFECT.contains(snpEffSummary.canonicalCodingEffect());
            }
        }

        return false;
    }

    @Override
    public void flush() {

    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template) {
        template.addMetaDataLine(new VCFInfoHeaderLine(REPORTED_FLAG, 0, VCFHeaderLineType.Flag, REPORTED_DESC));
        return template;
    }
}
