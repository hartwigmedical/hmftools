package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.apiclients.civic.api.CivicApiWrapper;
import com.hartwig.hmftools.apiclients.civic.data.CivicVariant;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.variant.Variant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.patientreporter.HmfReporterData;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;

@Value.Immutable
@Value.Style(allParameters = true)
public abstract class VariantReporterData {
    public static final FieldBuilder<?> VARIANT_NAME = field("variantName", String.class);
    public static final FieldBuilder<?> EVIDENCE_LINK = field("evidence_link", String.class);

    @NotNull
    public abstract String getVariantName();

    @NotNull
    public abstract List<CivicVariantReporterData> getVariants();

    public static List<VariantReporterData> of(@NotNull final PatientReport report, @NotNull final HmfReporterData reporterData) {
        final List<VariantReporterData> variantReporterData = Lists.newArrayList();
        for (final VariantReport variantReport : report.variants()) {
            for (final HmfGenomeRegion region : reporterData.slicer().hmfRegions()) {
                if (region.gene().equals(variantReport.gene())) {
                    final int entrezId = Integer.parseInt(region.entrezId());
                    final Variant variant = variantReportToVariant(variantReport);
                    final String variantName =
                            entrezId + "(" + variantReport.gene() + ")" + "\t" + variantReport.chromosomePosition() + "\t"
                                    + variantReport.variantField();
                    final List<CivicVariant> civicVariants = CivicApiWrapper.getVariantsContaining(entrezId, variant)
                            .filter(civicVariant -> !civicVariant.evidenceItemsWithDrugs().isEmpty())
                            .toList()
                            .blockingGet();
                    final List<CivicVariantReporterData> civicVariantReporterData =
                            civicVariants.stream().map(CivicVariantReporterData::of).collect(Collectors.toList());
                    if (!civicVariantReporterData.isEmpty()) {
                        variantReporterData.add(ImmutableVariantReporterData.of(variantName, civicVariantReporterData));
                    }
                }
            }
        }
        return variantReporterData;
    }

    public static Variant variantReportToVariant(@NotNull final VariantReport variantReport) {
        return new Variant() {
            @NotNull
            @Override
            public String ref() {
                return variantReport.ref();
            }

            @NotNull
            @Override
            public String alt() {
                return variantReport.alt();
            }

            @NotNull
            @Override
            public VariantType type() {
                return null;
            }

            @NotNull
            @Override
            public String filter() {
                return null;
            }

            @NotNull
            @Override
            public String chromosome() {
                return variantReport.chromosome();
            }

            @Override
            public long position() {
                return variantReport.position();
            }
        };
    }
}
