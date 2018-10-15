package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.cnv.ActionabilityCNVs;
import com.hartwig.hmftools.common.actionability.cnv.ActionabilityCNVsEvidenceItems;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.dynamicreports.report.definition.ReportParameters;
import net.sf.jasperreports.engine.JRDataSource;

public class EvidenceItemDataSourceCNVs {

    public static final FieldBuilder<?> GENE = field("Gene", String.class);
    public static final FieldBuilder<?> CHROMOSOME = field("Chromosome", String.class);
    public static final FieldBuilder<?> CHROMSOME_BAND = field("Chromome band", String.class);
    public static final FieldBuilder<?> CNV_TYPE= field("Cnv type", String.class);
    public static final FieldBuilder<?> DRUG = field("drug", String.class);
    public static final FieldBuilder<?> DRUGS_TYPE = field("drugs type", String.class);
    public static final FieldBuilder<?> LEVEL = field("level", String.class);
    public static final FieldBuilder<?> RESPONSE = field("response", String.class);
    public static final FieldBuilder<?> SOURCE = field("source", String.class);
    public static final FieldBuilder<?> LABEL = field("label", String.class);

    private EvidenceItemDataSourceCNVs() {
    }

    @NotNull
    public static JRDataSource fromActionabilityCNV(@NotNull Map<GeneCopyNumber, ActionabilityCNVsEvidenceItems> evidenceItemCNV) {
        final DRDataSource actionabilityCNVDatasource = new DRDataSource(GENE.getName(),
                CHROMOSOME.getName(),
                CHROMSOME_BAND.getName(),
                CNV_TYPE.getName(),
                DRUG.getName(),
                DRUGS_TYPE.getName(),
                LEVEL.getName(),
                RESPONSE.getName(),
                SOURCE.getName(),
                LABEL.getName());

        for (Map.Entry<GeneCopyNumber, ActionabilityCNVsEvidenceItems> entry : evidenceItemCNV.entrySet()) {

            String chromosome = entry.getKey().chromosome();
            String chromosomeBand = entry.getKey().chromosomeBand();

            for (ActionabilityCNVs CNV : entry.getValue().onLabel()) {
                actionabilityCNVDatasource.add(CNV.gene(),
                        chromosome,
                        chromosomeBand,
                        CNV.cnvType(),
                        CNV.drugsName(),
                        CNV.drugsType(),
                        CNV.hmfLevel(),
                        CNV.hmfResponse(),
                        sourceName(CNV.source()),
                        "yes");
            }

            for (ActionabilityCNVs CNV : entry.getValue().onLabel()) {
                actionabilityCNVDatasource.add(CNV.gene(),
                        chromosome,
                        chromosomeBand,
                        CNV.cnvType(),
                        CNV.drugsName(),
                        CNV.drugsType(),
                        CNV.hmfLevel(),
                        CNV.hmfResponse(),
                        sourceName(CNV.source()),
                        "no");
            }
        }
        return actionabilityCNVDatasource;
    }

    @NotNull
    private static String sourceName(@NotNull String source) {
        String sourceName = "";
        if (source.equals("oncoKb")) {
            sourceName = "OncoKB";
        } else if (source.equals("iclusion")) {
            sourceName = "Iclusion";
        } else if (source.equals("civic")) {
            sourceName = "CiViC";
        } else if (source.equals("cgi")) {
            sourceName = "CGI";
        }
        return sourceName;
    }

    @NotNull
    public static AbstractSimpleExpression<String> sourceHyperlink(@NotNull List<ActionabilityCNVs> evidenceItems) {
        List<String> linkReference = Lists.newArrayList();
        List<String> linkGene = Lists.newArrayList();
        for (ActionabilityCNVs variant : evidenceItems) {
            linkReference.add(variant.reference());
            linkGene.add(variant.gene());
        }
        return new AbstractSimpleExpression<String>() {
            @Override
            public String evaluate(@NotNull final ReportParameters data) {
                final String source = data.getValue(SOURCE.getName());
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
                        return "";
                }
            }
        };
    }

    @NotNull
    public static FieldBuilder<?>[] actionabilityCNVFields() {
        return new FieldBuilder<?>[] { GENE, CHROMOSOME, CHROMSOME_BAND, CNV_TYPE, DRUG, DRUGS_TYPE, LEVEL, RESPONSE, SOURCE, LABEL };
    }

}
