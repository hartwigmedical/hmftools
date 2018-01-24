package com.hartwig.hmftools.patientreporter.report.components;

import static com.hartwig.hmftools.patientreporter.report.Commons.HEADER_TO_DETAIL_VERTICAL_GAP;
import static com.hartwig.hmftools.patientreporter.report.Commons.baseTable;
import static com.hartwig.hmftools.patientreporter.report.Commons.dataStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.linkStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.toList;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.hyperLink;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientreporter.HmfReporterData;
import com.hartwig.hmftools.patientreporter.report.data.GenePanelDataSource;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.builder.component.VerticalListBuilder;

public final class GenePanelSection {

    @NotNull
    public static ComponentBuilder<?, ?> build(@NotNull final HmfReporterData reporterData) {
        final long coverage = Math.round(reporterData.panelGeneModel().numberOfBases() / 1E6);
        final VerticalListBuilder section = toList("Details on the reported gene panel",
                Lists.newArrayList("Findings are reported for the " + Integer.toString(reporterData.panelGeneModel().numberOfRegions())
                        + " genes (canonical transcripts) indicated below, covering " + coverage + " MBases."));

        return section.add(cmp.verticalGap(HEADER_TO_DETAIL_VERTICAL_GAP), genePanelTable(reporterData));
    }

    @NotNull
    private static ComponentBuilder<?, ?> genePanelTable(@NotNull final HmfReporterData reporterData) {
        // KODU: Overwrite default font size to make the panel fit on one page.
        final int fontSize = 7;
        return cmp.subreport(baseTable().setColumnStyle(dataStyle().setFontSize(fontSize))
                .setPageColumnsPerPage(2)
                .columns(col.emptyColumn().setFixedWidth(20),
                        col.column("Gene", GenePanelDataSource.GENE_FIELD).setWidth(50),
                        col.column("Transcript", GenePanelDataSource.TRANSCRIPT_FIELD)
                                .setHyperLink(hyperLink(GenePanelDataSource.transcriptUrl()))
                                .setStyle(linkStyle().setFontSize(fontSize)),
                        col.column("CosmicGenes Type", GenePanelDataSource.TYPE_FIELD),
                        col.emptyColumn().setFixedWidth(20))).setDataSource(GenePanelDataSource.fromHmfReporterData(reporterData));
    }
}
