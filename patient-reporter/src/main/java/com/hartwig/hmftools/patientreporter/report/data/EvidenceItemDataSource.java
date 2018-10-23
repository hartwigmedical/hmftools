package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.actionability.EvidenceItem;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.dynamicreports.report.definition.ReportParameters;
import net.sf.jasperreports.engine.JRDataSource;

public final class EvidenceItemDataSource {

    public static final FieldBuilder<?> EVENT_FIELD = field("event", String.class);
    public static final FieldBuilder<?> DRUG_FIELD = field("drug", String.class);
    public static final FieldBuilder<?> DRUGS_TYPE_FIELD = field("drug type", String.class);
    public static final FieldBuilder<?> LEVEL_FIELD = field("level", String.class);
    public static final FieldBuilder<?> RESPONSE_FIELD = field("response", String.class);
    public static final FieldBuilder<?> SOURCE_FIELD = field("source", String.class);
    public static final FieldBuilder<?> ON_LABEL_FIELD = field("on_label", String.class);
    private static final FieldBuilder<?> REFERENCE_FIELD = field("reference", String.class);

    private EvidenceItemDataSource() {
    }

    @NotNull
    public static FieldBuilder<?>[] evidenceItemFields() {
        return new FieldBuilder<?>[] { EVENT_FIELD, DRUG_FIELD, DRUGS_TYPE_FIELD, LEVEL_FIELD, RESPONSE_FIELD, SOURCE_FIELD, ON_LABEL_FIELD,
                REFERENCE_FIELD };
    }

    @NotNull
    public static JRDataSource fromEvidenceItems(@NotNull List<EvidenceItem> evidenceItems) {
        final DRDataSource evidenceItemDataSource = new DRDataSource(EVENT_FIELD.getName(),
                DRUG_FIELD.getName(),
                DRUGS_TYPE_FIELD.getName(),
                LEVEL_FIELD.getName(),
                RESPONSE_FIELD.getName(),
                SOURCE_FIELD.getName(),
                REFERENCE_FIELD.getName(),
                ON_LABEL_FIELD.getName());

        for (EvidenceItem evidenceItem : sort(evidenceItems)) {
            if(!evidenceItem.source().isTrialSource() ){
                evidenceItemDataSource.add(evidenceItem.event(),
                        evidenceItem.drug(),
                        evidenceItem.drugsType(),
                        evidenceItem.level(),
                        evidenceItem.response(),
                        evidenceItem.source().sourceName(),
                        evidenceItem.reference(),
                        evidenceItem.isOnLabel() ? "Yes" : "No");
            }
        }
        return evidenceItemDataSource;
    }

    @NotNull
    private static List<EvidenceItem> sort(@NotNull List<EvidenceItem> evidenceItems) {
        return evidenceItems.stream().sorted((item1, item2) -> {
            if (item1.event().equals(item2.event())) {
                return item1.event().compareTo(item2.event());
            } else {
                return item1.level().compareTo(item2.level());
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    public static AbstractSimpleExpression<String> sourceHyperlink() {
        return new AbstractSimpleExpression<String>() {
            @Override
            public String evaluate(@NotNull final ReportParameters data) {
                String source = data.getValue(SOURCE_FIELD.getName()).toString();
                String reference = data.getValue(REFERENCE_FIELD.getName()).toString();
                String gene = data.getValue(EVENT_FIELD.getName()).toString();
                switch (source.toLowerCase()) {
                    case "oncokb":
                        String[] geneId = gene.split(" ");
                        return "http://oncokb.org/#/gene/" + geneId[0] + "/alteration/" + reference;
                    case "cgi":
                        return "https://www.cancergenomeinterpreter.org/biomarkers";
                    case "civic":
                        String[] variantId = reference.split(":");
                        return "https://civic.genome.wustl.edu/links/variants/" + variantId[1];
                    default:
                        return Strings.EMPTY;
                }
            }
        };
    }
}
