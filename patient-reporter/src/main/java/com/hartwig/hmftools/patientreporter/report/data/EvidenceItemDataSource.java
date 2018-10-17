package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;

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

public abstract class EvidenceItemDataSource {
    private static final Logger LOGGER = LogManager.getLogger(EvidenceItemDataSource.class);

    public static final FieldBuilder<?> GENE_FIELD = field("Gene", String.class);
    public static final FieldBuilder<?> EVENT_FIELD = field("variant", String.class);
    public static final FieldBuilder<?> DRUG_FIELD = field("drug", String.class);
    public static final FieldBuilder<?> DRUGS_TYPE_FIELD = field("drug type", String.class);
    public static final FieldBuilder<?> LEVEL_FIELD = field("level", String.class);
    public static final FieldBuilder<?> RESPONSE_FIELD = field("response", String.class);
    public static final FieldBuilder<?> SOURCE_FIELD = field("source", String.class);
    private static final FieldBuilder<?> REFERENCE_FIELD = field("reference", String.class);
    public static final FieldBuilder<?> ON_LABEL_FIELD = field("on_label", String.class);

    private EvidenceItemDataSource() {
    }

    @NotNull
    public static FieldBuilder<?>[] evidenceItemFields() {
        return new FieldBuilder<?>[] { GENE_FIELD, EVENT_FIELD, DRUG_FIELD, DRUGS_TYPE_FIELD, LEVEL_FIELD, RESPONSE_FIELD, SOURCE_FIELD,
                REFERENCE_FIELD, ON_LABEL_FIELD };
    }

    @NotNull
    public static JRDataSource fromEvidenceItems(@NotNull List<EvidenceItem> evidenceItems) {
        final DRDataSource evidenceItemDataSource = new DRDataSource(GENE_FIELD.getName(),
                EVENT_FIELD.getName(),
                DRUG_FIELD.getName(),
                DRUGS_TYPE_FIELD.getName(),
                LEVEL_FIELD.getName(),
                RESPONSE_FIELD.getName(),
                SOURCE_FIELD.getName(),
                REFERENCE_FIELD.getName(),
                ON_LABEL_FIELD.getName());

        for (EvidenceItem evidenceItem : evidenceItems) {
            evidenceItemDataSource.add(evidenceItem.gene(),
                    evidenceItem.event(),
                    evidenceItem.drug(),
                    evidenceItem.drugsType(),
                    evidenceItem.level(),
                    evidenceItem.response(),
                    sourceName(evidenceItem.source()),
                    evidenceItem.reference(),
                    evidenceItem.isOnLabel() ? "Yes" : "No");
        }

        return evidenceItemDataSource;
    }

    @NotNull
    private static String sourceName(@NotNull String source) {
        switch (source) {
            case "oncoKb":
                return "OncoKB";
            case "iclusion":
                return "Iclusion";
            case "civic":
                return "CiViC";
            case "cgi":
                return "CGI";
            default:
                LOGGER.warn("Unrecognized source in evidence item: " + source);
                return Strings.EMPTY;
        }
    }

    @NotNull
    public static AbstractSimpleExpression<String> sourceHyperlink() {
        return new AbstractSimpleExpression<String>() {
            @Override
            public String evaluate(@NotNull final ReportParameters data) {
                String source = data.getValue(SOURCE_FIELD.getName()).toString().toLowerCase();
                // TODO (LISC): Use reference to create link.
                String reference = data.getValue(REFERENCE_FIELD.getName()).toString();
                switch (source) {
                    case "oncoKb":
                        return "http://oncokb.org/#/gene/";
                    //+ linkGene.get(0) + "/alteration/" + linkReference.get(0);
                    case "iclusion":
                        return "https://www.iclusion.org";
                    case "cgi":
                        return "https://www.cancergenomeinterpreter.org/biomarkers";
                    case "civic":
                        //  String[] link = linkReference.get(0).split(":");
                        return "https://civic.genome.wustl.edu/links/variants/";
                    default:
                        return Strings.EMPTY;
                }
            }
        };
    }
}
