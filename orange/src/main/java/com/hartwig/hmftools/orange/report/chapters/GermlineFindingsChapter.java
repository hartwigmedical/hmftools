package com.hartwig.hmftools.orange.report.chapters;

import java.text.DecimalFormat;
import java.util.Map;

import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.tables.GermlineVariantTable;
import com.hartwig.hmftools.orange.report.util.TableUtil;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;

public class GermlineFindingsChapter implements ReportChapter {

    private static final DecimalFormat PERCENTAGE_FORMAT = ReportResources.decimalFormat("#.0'%'");

    @NotNull
    private final OrangeReport report;
    private final boolean reportGermline;

    public GermlineFindingsChapter(@NotNull final OrangeReport report, final boolean reportGermline) {
        this.report = report;
        this.reportGermline = reportGermline;
    }

    @NotNull
    @Override
    public String name() {
        return "Germline Findings";
    }

    @Override
    public void render(@NotNull final Document document) {
        document.add(new Paragraph("Germline Findings"));
        if (reportGermline) {
            addGermlineVariants(document);
            addMVLHAnalysis(document);
            addGermlineCNAberrations(document);
            addPEACH(document);
        } else {
            document.add(new Paragraph("NA").addStyle(ReportResources.tableContentStyle()));
        }
    }

    private void addGermlineVariants(@NotNull Document document) {
        String germlineDriversTitle = "Driver variants (" + report.purple().reportableGermlineVariants().size() + ")";
        Table germlineDriversTable = GermlineVariantTable.build(germlineDriversTitle, report.purple().reportableGermlineVariants());
        document.add(germlineDriversTable);

        document.add(new Paragraph("TODO: Add other potentially other germline variants").addStyle(ReportResources.tableContentStyle()));
    }

    private void addMVLHAnalysis(@NotNull Document document) {
        Table table = TableUtil.createReportContentTable(new float[] { 2, 2, 1, 2, 2, 1, 2, 2, 1, 2, 2, 1 },
                new Cell[] { TableUtil.createHeaderCell("Gene"), TableUtil.createHeaderCell("MVLH"), TableUtil.createHeaderCell(""),
                        TableUtil.createHeaderCell("Gene"), TableUtil.createHeaderCell("MVLH"), TableUtil.createHeaderCell(""),
                        TableUtil.createHeaderCell("Gene"), TableUtil.createHeaderCell("MVLH"), TableUtil.createHeaderCell(""),
                        TableUtil.createHeaderCell("Gene"), TableUtil.createHeaderCell("MVLH"), TableUtil.createHeaderCell("")});

        for (Map.Entry<String, Double> mvlhEntry : report.germlineMVLHPerGene().entrySet()) {
            double mvlh = mvlhEntry.getValue();
            if (mvlh > 0.01) {
                table.addCell(TableUtil.createContentCell(mvlhEntry.getKey()));
                table.addCell(TableUtil.createContentCell(PERCENTAGE_FORMAT.format(mvlh * 100)));
                table.addCell(TableUtil.createContentCell(""));
            }
        }
        document.add(TableUtil.createWrappingReportTable(table, "Genes with missed variant likelihood > 1%"));
    }

    private void addGermlineCNAberrations(@NotNull Document document) {
        document.add(new Paragraph("TODO: Add Germline CN aberrations").addStyle(ReportResources.tableContentStyle()));
    }

    private void addPEACH(@NotNull Document document) {
        document.add(new Paragraph("TODO: Add PEACH").addStyle(ReportResources.tableContentStyle()));
    }

}
