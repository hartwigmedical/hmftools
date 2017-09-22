package com.hartwig.hmftools.patientreporter.report.data;

import static com.hartwig.hmftools.patientreporter.report.data.EvidenceReportData.variantReportToVariant;

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
public abstract class Alteration {
    private static final Logger LOGGER = LogManager.getLogger(Alteration.class);

    public static final FieldBuilder<?> ALTERATION = field("alteration", String.class);

    public abstract String getGene();

    public abstract String getPredictedEffect();

    public abstract List<AlterationEvidence> getEvidence();

    public abstract List<AlterationMatch> getMatches();

    public String getAlteration() {
        return getGene() + "\n" + getPredictedEffect();
    }

    public static Alteration from(@NotNull final VariantReport variantReport, @NotNull final List<CivicVariant> civicVariants) {
        final String gene = variantReport.gene();
        final String predictedEffect = variantReport.hgvsProtein();
        final List<AlterationEvidence> exactMatchEvidence = Lists.newArrayList();
        final List<AlterationMatch> matchingVariants = Lists.newArrayList();

        civicVariants.forEach(civicVariant -> {
            if (civicVariant.coordinates().equals(variantReportToVariant(variantReport))) {
                for (final String significance : civicVariant.groupedEvidenceItems().keySet()) {
                    final String drugsString = getDrugsWithEvidenceLevel(civicVariant.groupedEvidenceItems().get(significance));
                    exactMatchEvidence.add(ImmutableAlterationEvidence.of(significance, drugsString, "CIViC"));
                }
                matchingVariants.add(AlterationMatch.of("exact", civicVariant));
            } else {
                matchingVariants.add(AlterationMatch.of("approx.", civicVariant));
            }
        });
        return ImmutableAlteration.of(gene, predictedEffect, exactMatchEvidence, matchingVariants);
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
