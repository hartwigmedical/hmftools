package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.apiclients.civic.data.CivicEvidenceItem;
import com.hartwig.hmftools.apiclients.civic.data.CivicVariant;

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
    public static final FieldBuilder<?> LOH = field("loh", String.class);
    public static final FieldBuilder<?> SUBCLONAL = field("subclonal", String.class);

    public abstract String getGene();

    public abstract String getPredictedEffect();

    public abstract String getLoh();

    public abstract String getSubclonal();

    public abstract List<AlterationEvidenceReporterData> getEvidence();

    public String getAlteration() {
        return getGene() + "\n" + getPredictedEffect();
    }

    public static AlterationReporterData from(@NotNull final String gene, @NotNull final String predictedEffect, @NotNull final String loh,
            @NotNull final String subclonal, @NotNull final List<CivicVariant> civicVariants) {
        if (civicVariants.size() > 1) {
            LOGGER.warn("found more than one civic variant for: " + gene + "(" + predictedEffect + ")");
        }
        final List<AlterationEvidenceReporterData> alterationEvidence = Lists.newArrayList();
        civicVariants.forEach(civicVariant -> {
            for (final String significance : civicVariant.groupedEvidenceItems().keySet()) {
                final List<String> drugs = civicVariant.groupedEvidenceItems()
                        .get(significance)
                        .keySet()
                        .stream()
                        .map(drug -> drug + "(" + Strings.join(civicVariant.groupedEvidenceItems()
                                .get(significance)
                                .get(drug)
                                .stream()
                                .map(CivicEvidenceItem::level)
                                .distinct()
                                .collect(Collectors.toList()), ',') + ")")
                        .collect(Collectors.toList());
                final String drugsString = Strings.join(drugs, '\n');
                alterationEvidence.add(ImmutableAlterationEvidenceReporterData.of(significance, drugsString, "CIViC"));
            }
        });
        return ImmutableAlterationReporterData.of(gene, predictedEffect, loh, subclonal, alterationEvidence);
    }
}
