package com.hartwig.hmftools.patientreporter.report.pages;

import static com.hartwig.hmftools.patientreporter.report.Commons.SECTION_VERTICAL_GAP;
import static com.hartwig.hmftools.patientreporter.report.Commons.fontStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.sectionHeaderStyle;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.exp;
import static net.sf.dynamicreports.report.builder.DynamicReports.hyperLink;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;

import java.awt.Color;

import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.report.Commons;
import com.hartwig.hmftools.patientreporter.report.components.MainPageTopSection;
import com.hartwig.hmftools.patientreporter.report.data.Alteration;
import com.hartwig.hmftools.patientreporter.report.data.AlterationEvidence;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.builder.component.SubreportBuilder;
import net.sf.dynamicreports.report.builder.style.StyleBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;

@Value.Immutable
@Value.Style(passAnnotations = NotNull.class,
             allParameters = true)
public abstract class AlterationEvidencePage {
    private static final int PADDING = 3;

    @NotNull
    public static ComponentBuilder<?, ?> reportComponent(@NotNull final SampleReport sampleReport) {
        return cmp.verticalList(MainPageTopSection.build("HMF Civic Evidence Supplement", sampleReport),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("Knowledgebase drug association of reported genomic alterations").setStyle(sectionHeaderStyle().setFontSize(15)),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                conciseEvidenceSection());
    }

    @NotNull
    private static ComponentBuilder<?, ?> alterationEvidenceTable() {
        final int ALTERATION_WIDTH = 130;
        final int SIGNIFICANCE_WIDTH = 100;
        final int SOURCE_WIDTH = 70;

        final SubreportBuilder subTable = cmp.subreport(alterationBaseTable().setColumnStyle(alterationDataStyle())
                .fields(AlterationEvidence.SOURCE_URL)
                .columns(col.column(AlterationEvidence.SIGNIFICANCE).setFixedWidth(SIGNIFICANCE_WIDTH).setMinHeight(25),
                        col.column(AlterationEvidence.DRUGS),
                        col.column(AlterationEvidence.SOURCE)
                                .setHyperLink(hyperLink(AlterationEvidence.sourceHyperlink()))
                                .setStyle(alterationDataLinkStyle())
                                .setFixedWidth(SOURCE_WIDTH))).setDataSource(exp.subDatasourceBeanCollection("evidence"));

        final ComponentBuilder<?, ?> tableHeader =
                cmp.horizontalList(cmp.text("Alteration").setStyle(alterationTableHeaderStyle()).setFixedWidth(ALTERATION_WIDTH),
                        cmp.text("Significance").setStyle(alterationTableHeaderStyle()).setFixedWidth(SIGNIFICANCE_WIDTH),
                        cmp.text("Association(Lv)").setStyle(alterationTableHeaderStyle()),
                        cmp.text("Source").setStyle(alterationTableHeaderStyle()).setFixedWidth(SOURCE_WIDTH));

        return cmp.subreport(alterationBaseTable().setColumnStyle(alterationDataStyle())
                .title(tableHeader)
                .columns(col.column(Alteration.ALTERATION).setFixedWidth(ALTERATION_WIDTH), col.componentColumn(subTable))
                .noData(cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER))))
                .setDataSource(exp.subDatasourceBeanCollection("alterationsWithEvidence"));
    }

    @NotNull
    private static ComponentBuilder<?, ?> conciseEvidenceSection() {
        return cmp.horizontalList(cmp.horizontalGap(20), alterationEvidenceTable(), cmp.horizontalGap(20));
    }

    @NotNull
    static JasperReportBuilder alterationBaseTable() {
        return report().setColumnStyle(alterationDataStyle()).setColumnTitleStyle(alterationTableHeaderStyle());
    }

    @NotNull
    static StyleBuilder alterationTableHeaderStyle() {
        return Commons.tableHeaderStyle().setFontSize(9).setPadding(PADDING);
    }

    @NotNull
    static StyleBuilder alterationDataStyle() {
        return Commons.smallDataTableStyle().setPadding(PADDING);
    }

    @NotNull
    static StyleBuilder alterationDataLinkStyle() {
        return alterationDataStyle().setForegroundColor(Color.BLUE);
    }
}
