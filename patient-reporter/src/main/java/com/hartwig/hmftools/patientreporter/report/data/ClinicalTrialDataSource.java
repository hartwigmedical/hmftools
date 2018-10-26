package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.actionability.EvidenceItem;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.dynamicreports.report.definition.ReportParameters;
import net.sf.jasperreports.engine.JRDataSource;

public class ClinicalTrialDataSource {

    public static final FieldBuilder<?> EVENT_FIELD = field("event", String.class);
    public static final FieldBuilder<?> TRIAL_FIELD = field("trial", String.class);
    public static final FieldBuilder<?> SOURCE_FIELD = field("source", String.class);
    public static final FieldBuilder<?> CCMO_FIELD = field("ccmo", String.class);
    public static final FieldBuilder<?> ON_LABEL_FIELD = field("on_label", String.class);
    private static final FieldBuilder<?> REFERENCE_FIELD = field("reference", String.class);

    private ClinicalTrialDataSource() {
    }

    @NotNull
    public static FieldBuilder<?>[] clinicalTrialFields() {
        return new FieldBuilder<?>[] { EVENT_FIELD, TRIAL_FIELD, SOURCE_FIELD, CCMO_FIELD, ON_LABEL_FIELD, REFERENCE_FIELD };
    }

    @NotNull
    public static JRDataSource fromClinicalTrials(@NotNull List<EvidenceItem> trials) {
        final DRDataSource evidenceItemDataSource = new DRDataSource(EVENT_FIELD.getName(),
                TRIAL_FIELD.getName(),
                SOURCE_FIELD.getName(),
                CCMO_FIELD.getName(),
                ON_LABEL_FIELD.getName(),
                REFERENCE_FIELD.getName());

        for (EvidenceItem evidenceItem : sort(trials)) {
            assert evidenceItem.source().isTrialSource();

            evidenceItemDataSource.add(evidenceItem.event(),
                    evidenceItem.drug(),
                    evidenceItem.source().sourceName(),
                    CCMOId(evidenceItem.reference()),
                    evidenceItem.isOnLabel() ? "Yes" : "No",
                    evidenceItem.reference());
        }

        return evidenceItemDataSource;
    }

    @NotNull
    private static String CCMOId(@NotNull String reference) {
        // KODU: Expected format "EXT1 (CCMO)"
        String referenceWithoutParenthesis = reference.replace(")", "");
        String[] splitExtAndCCMO = referenceWithoutParenthesis.split("\\(");
        return splitExtAndCCMO[1];
    }

    @NotNull
    private static List<EvidenceItem> sort(@NotNull List<EvidenceItem> evidenceItems) {
        return evidenceItems.stream().sorted(Comparator.comparing(EvidenceItem::drug)).collect(Collectors.toList());
    }

    @NotNull
    public static AbstractSimpleExpression<String> sourceHyperlink() {
        return new AbstractSimpleExpression<String>() {
            @Override
            public String evaluate(@NotNull final ReportParameters data) {
                String source = data.getValue(SOURCE_FIELD.getName()).toString();
                String reference = data.getValue(REFERENCE_FIELD.getName()).toString();
                String ext = EXTId(reference);
                switch (source.toLowerCase()) {
                    case "iclusion":
                        return "https://iclusion.org/hmf/" + ext;
                    default:
                        return Strings.EMPTY;
                }
            }
        };
    }

    @NotNull
    private static String EXTId(@NotNull String reference) {
        // KODU: Expected format "EXT1 (CCMO)"
        String[] splitExtAndCCMO = reference.split("\\(");
        String ext = splitExtAndCCMO[0];
        return ext.substring(3, ext.length()).trim();
    }
}
