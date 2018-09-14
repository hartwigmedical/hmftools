package com.hartwig.hmftools.patientreporter.report.pages;

import static com.hartwig.hmftools.patientreporter.report.Commons.HEADER_TO_DETAIL_VERTICAL_GAP;
import static com.hartwig.hmftools.patientreporter.report.Commons.SECTION_VERTICAL_GAP;
import static com.hartwig.hmftools.patientreporter.report.Commons.baseTable;
import static com.hartwig.hmftools.patientreporter.report.Commons.dataStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.linkStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.sectionHeaderStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.toList;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.hyperLink;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientreporter.SequencedReportData;
import com.hartwig.hmftools.patientreporter.report.Commons;
import com.hartwig.hmftools.patientreporter.report.data.GenePanelDataSource;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.builder.component.VerticalListBuilder;

@Value.Immutable
@Value.Style(passAnnotations = NotNull.class,
             allParameters = true)
public abstract class GenePanelPage {

    @NotNull
    abstract SequencedReportData reporterData();

    public ComponentBuilder<?, ?> reportComponent() {
        return cmp.verticalList(cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text(Commons.TITLE_SEQUENCE+ " - Gene Panel").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                build(reporterData()));
    }

    @NotNull
    private static ComponentBuilder<?, ?> build(@NotNull final SequencedReportData reporterData) {
        final long coverage = Math.round(reporterData.panelGeneModel().numberOfBases() / 1E6);
        final VerticalListBuilder section = toList("Details on the reported gene panel",
                Lists.newArrayList("Findings are reported for the " + Integer.toString(reporterData.panelGeneModel().numberOfRegions())
                        + " genes (canonical transcripts) indicated below, covering " + coverage + " MBases."));

        return section.add(cmp.verticalGap(HEADER_TO_DETAIL_VERTICAL_GAP), genePanelTable(reporterData));
    }

    @NotNull
    private static ComponentBuilder<?, ?> genePanelTable(@NotNull final SequencedReportData reporterData) {
        // KODU: Overwrite default font size to make the panel fit on one page.
        final int fontSize = 6;
        return cmp.subreport(baseTable().setColumnStyle(dataStyle().setFontSize(fontSize))
                .setPageColumnsPerPage(2)
                .columns(col.emptyColumn().setFixedWidth(20),
                        col.column("Gene", GenePanelDataSource.GENE_FIELD).setWidth(50),
                        col.column("Transcript", GenePanelDataSource.TRANSCRIPT_FIELD)
                                .setHyperLink(hyperLink(GenePanelDataSource.transcriptUrl()))
                                .setStyle(linkStyle().setFontSize(fontSize)),
                        col.column("Type", GenePanelDataSource.TYPE_FIELD),
                        col.emptyColumn().setFixedWidth(20))).setDataSource(GenePanelDataSource.fromSequencedReportData(reporterData));
    }
}
