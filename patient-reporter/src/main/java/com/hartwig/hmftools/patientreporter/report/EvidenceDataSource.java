package com.hartwig.hmftools.patientreporter.report;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.civic.data.CivicDrug;
import com.hartwig.hmftools.civic.data.CivicEvidenceItem;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

class EvidenceDataSource {
    static final FieldBuilder<?> LEVEL_FIELD = field("level", Character.class);
    static final FieldBuilder<?> TUMOR_TYPE_FIELD = field("tumorType", String.class);
    static final FieldBuilder<?> DIRECTION_FIELD = field("direction", String.class);
    static final FieldBuilder<?> SIGNIFICANCE_FIELD = field("significance", String.class);
    static final FieldBuilder<?> DRUGS_FIELD = field("drugs", String.class);

    private EvidenceDataSource() {
    }

    @NotNull
    static JRDataSource fromEvidenceItems(@NotNull final List<CivicEvidenceItem> evidenceItems) {
        final DRDataSource genePanelDataSource =
                new DRDataSource(LEVEL_FIELD.getName(), TUMOR_TYPE_FIELD.getName(), DIRECTION_FIELD.getName(), SIGNIFICANCE_FIELD.getName(),
                        DRUGS_FIELD.getName());
        for (final CivicEvidenceItem evidenceItem : evidenceItems) {
            final String drugs = Strings.join(evidenceItem.drugs().stream().map(CivicDrug::name).collect(Collectors.toList()), ',');

            genePanelDataSource.add(evidenceItem.level(), evidenceItem.disease().displayName(), evidenceItem.direction(),
                    evidenceItem.significance(), drugs);
        }
        return genePanelDataSource;
    }

    @NotNull
    static JRDataSource fromGroupedEvidenceItems(@NotNull final Map<String, Map<String, List<CivicEvidenceItem>>> groupedEvidenceItems) {
        final DRDataSource genePanelDataSource = new DRDataSource(SIGNIFICANCE_FIELD.getName(), DRUGS_FIELD.getName());
        for (final String significance : groupedEvidenceItems.keySet()) {
            final List<String> drugs = groupedEvidenceItems.get(significance)
                    .keySet()
                    .stream()
                    .map(drug -> drug + "(" + Strings.join(groupedEvidenceItems.get(significance)
                            .get(drug)
                            .stream()
                            .map(CivicEvidenceItem::level)
                            .distinct()
                            .collect(Collectors.toList()), ',') + ")")
                    .collect(Collectors.toList());
            final String drugsString = Strings.join(drugs, ',');
            genePanelDataSource.add(significance, drugsString);
        }
        return genePanelDataSource;
    }

    @NotNull
    static FieldBuilder<?>[] evidenceFields() {
        return new FieldBuilder<?>[] { LEVEL_FIELD, TUMOR_TYPE_FIELD, DIRECTION_FIELD, SIGNIFICANCE_FIELD, DRUGS_FIELD };
    }

    @NotNull
    static FieldBuilder<?>[] conciseEvidenceFields() {
        return new FieldBuilder<?>[] { SIGNIFICANCE_FIELD, DRUGS_FIELD };
    }
}
