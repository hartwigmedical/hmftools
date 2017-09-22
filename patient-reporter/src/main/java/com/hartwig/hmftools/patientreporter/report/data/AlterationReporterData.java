package com.hartwig.hmftools.patientreporter.report.data;

import static com.hartwig.hmftools.patientreporter.report.data.VariantReporterData.variantReportToVariant;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.apiclients.civic.data.CivicEvidenceItem;
import com.hartwig.hmftools.apiclients.civic.data.CivicVariant;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;

@Value.Immutable
@Value.Style(allParameters = true)
public abstract class AlterationReporterData {
    private static final Logger LOGGER = LogManager.getLogger(AlterationReporterData.class);

    public static final FieldBuilder<?> ALTERATION = field("alteration", String.class);

    public abstract String getGene();

    public abstract String getPredictedEffect();

    public abstract List<AlterationEvidenceReporterData> getEvidence();

    public abstract List<AlterationMatch> getMatches();

    public String getAlteration() {
        return getGene() + "\n" + getPredictedEffect();
    }

    public static AlterationReporterData from(@NotNull final VariantReport variantReport, @NotNull final List<CivicVariant> civicVariants) {
        final String gene = variantReport.gene();
        final String predictedEffect = variantReport.hgvsProtein();
        final List<AlterationEvidenceReporterData> exactMatchEvidence = Lists.newArrayList();
        final List<AlterationMatch> matchingVariants = Lists.newArrayList();

        civicVariants.forEach(civicVariant -> {
            if (civicVariant.coordinates().equals(variantReportToVariant(variantReport))) {
                for (final String significance : civicVariant.groupedEvidenceItems().keySet()) {
                    final String drugsString = getDrugsWithEvidenceLevel(civicVariant.groupedEvidenceItems().get(significance));
                    exactMatchEvidence.add(ImmutableAlterationEvidenceReporterData.of(significance, drugsString, "CIViC"));
                }
                matchingVariants.add(AlterationMatch.of("exact", civicVariant));
            } else {
                matchingVariants.add(AlterationMatch.of("approx.", civicVariant));
            }
        });
        return ImmutableAlterationReporterData.of(gene, predictedEffect, exactMatchEvidence, matchingVariants);
    }

    @NotNull
    private static String getDrugsWithEvidenceLevel(@NotNull final Map<String, List<CivicEvidenceItem>> drugToEvidenceItems) {
        final List<String> drugs = drugToEvidenceItems.entrySet()
                .stream()
                .map(entry -> entry.getKey() + "(" + Strings.join(
                        entry.getValue().stream().map(CivicEvidenceItem::level).distinct().collect(Collectors.toList()), ',') + ")")
                .collect(Collectors.toList());
        return Strings.join(drugs, '\n');
    }
}
