package com.hartwig.hmftools.patientreporter.report.pages;

import static com.hartwig.hmftools.patientreporter.report.Commons.SECTION_VERTICAL_GAP;
import static com.hartwig.hmftools.patientreporter.report.Commons.fontStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.sectionHeaderStyle;
import static com.hartwig.hmftools.patientreporter.report.pages.AlterationEvidencePage.alterationBaseTable;
import static com.hartwig.hmftools.patientreporter.report.pages.AlterationEvidencePage.alterationDataLinkStyle;
import static com.hartwig.hmftools.patientreporter.report.pages.AlterationEvidencePage.alterationDataStyle;
import static com.hartwig.hmftools.patientreporter.report.pages.AlterationEvidencePage.alterationTableHeaderStyle;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.exp;
import static net.sf.dynamicreports.report.builder.DynamicReports.hyperLink;

import com.hartwig.hmftools.patientreporter.report.data.Alteration;
import com.hartwig.hmftools.patientreporter.report.data.AlterationMatch;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.builder.component.SubreportBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;

@Value.Immutable
@Value.Style(passAnnotations = NotNull.class,
             allParameters = true)
public abstract class AlterationDebugPage {

    @NotNull
    public static ComponentBuilder<?, ?> reportComponent() {
        return cmp.verticalList(cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("Civic matching variants").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                civicMatchingVariantsTable(),
                cmp.verticalGap(SECTION_VERTICAL_GAP));
    }

    @NotNull
    private static ComponentBuilder<?, ?> civicMatchingVariantsTable() {
        final int ALTERATION_WIDTH = 75;

        final SubreportBuilder subTable = cmp.subreport(alterationBaseTable().setColumnStyle(alterationDataStyle())
                .fields(AlterationMatch.SUMMARY_URL)
                .columns(col.column(AlterationMatch.MATCH_TYPE).setFixedWidth(40).setMinHeight(25),
                        col.column(AlterationMatch.NAME)
                                .setHyperLink(hyperLink(AlterationMatch.civicSummaryHyperlink()))
                                .setStyle(alterationDataLinkStyle())
                                .setFixedWidth(80),
                        col.column(AlterationMatch.VARIANT_TYPE).setFixedWidth(80),
                        col.column(AlterationMatch.CHROMOSOME).setFixedWidth(30),
                        col.column(AlterationMatch.START).setFixedWidth(45),
                        col.column(AlterationMatch.STOP).setFixedWidth(45),
                        col.column(AlterationMatch.CHROMOSOME2).setFixedWidth(30),
                        col.column(AlterationMatch.START2).setFixedWidth(30),
                        col.column(AlterationMatch.STOP2).setFixedWidth(30),
                        col.column(AlterationMatch.HGVS_EXPRESSIONS))).setDataSource(exp.subDatasourceBeanCollection("matches"));

        final ComponentBuilder<?, ?> tableHeader =
                cmp.horizontalList(cmp.text("Alteration").setStyle(alterationTableHeaderStyle()).setFixedWidth(ALTERATION_WIDTH),
                        cmp.text("Match").setStyle(alterationTableHeaderStyle()).setFixedWidth(40),
                        cmp.text("Name").setStyle(alterationTableHeaderStyle()).setFixedWidth(80),
                        cmp.text("Type").setStyle(alterationTableHeaderStyle()).setFixedWidth(80),
                        cmp.text("Chr").setStyle(alterationTableHeaderStyle()).setFixedWidth(30),
                        cmp.text("Start").setStyle(alterationTableHeaderStyle()).setFixedWidth(45),
                        cmp.text("Stop").setStyle(alterationTableHeaderStyle()).setFixedWidth(45),
                        cmp.text("Chr2").setStyle(alterationTableHeaderStyle()).setFixedWidth(30),
                        cmp.text("Start2").setStyle(alterationTableHeaderStyle()).setFixedWidth(30),
                        cmp.text("Stop2").setStyle(alterationTableHeaderStyle()).setFixedWidth(30),
                        cmp.text("Hgvs Expressions").setStyle(alterationTableHeaderStyle()));

        return cmp.subreport(alterationBaseTable().setColumnStyle(alterationDataStyle())
                .title(tableHeader)
                .columns(col.column(Alteration.ALTERATION).setFixedWidth(ALTERATION_WIDTH), col.componentColumn(subTable))
                .noData(cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER))))
                .setDataSource(exp.subDatasourceBeanCollection("alterations"));
    }
}
