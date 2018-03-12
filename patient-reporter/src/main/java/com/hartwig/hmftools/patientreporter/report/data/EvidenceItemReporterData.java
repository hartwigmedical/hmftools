package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.stream.Collectors;

import com.hartwig.hmftools.apiclients.civic.data.CivicDrug;
import com.hartwig.hmftools.apiclients.civic.data.CivicEvidenceItem;

import org.apache.logging.log4j.util.Strings;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;

@Value.Immutable
@Value.Style(allParameters = true)
public abstract class EvidenceItemReporterData {
    public static final FieldBuilder<?> LEVEL_FIELD = field("level", Character.class);
    public static final FieldBuilder<?> TUMOR_TYPE_FIELD = field("cancerType", String.class);
    public static final FieldBuilder<?> DIRECTION_FIELD = field("direction", String.class);
    public static final FieldBuilder<String> SIGNIFICANCE_FIELD = field("significance", String.class);
    public static final FieldBuilder<?> DRUGS_FIELD = field("drugs", String.class);

    public abstract Character getLevel();

    public abstract String getTumorType();

    public abstract String getDirection();

    public abstract String getSignificance();

    public abstract String getDrugs();

    public static EvidenceItemReporterData of(@NotNull final CivicEvidenceItem evidenceItem) {
        final String drugs = Strings.join(evidenceItem.drugs().stream().map(CivicDrug::name).collect(Collectors.toList()), ',');

        final String direction = evidenceItem.direction();
        final String directionString = direction == null ? "" : direction;

        final String significance = evidenceItem.significance();
        final String significanceString = significance == null ? "" : significance;

        return ImmutableEvidenceItemReporterData.of(evidenceItem.level(), evidenceItem.disease().displayName(), directionString,
                significanceString, drugs);
    }
}
