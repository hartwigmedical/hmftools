package com.hartwig.hmftools.orange.report.chapters;

import java.text.DecimalFormat;
import java.util.Comparator;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.tables.GermlineVariantTable;
import com.hartwig.hmftools.orange.report.util.TableUtil;
import com.itextpdf.kernel.geom.PageSize;
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

    @NotNull
    @Override
    public PageSize pageSize() {
        return PageSize.A4;
    }

    @Override
    public void render(@NotNull final Document document) {
        document.add(new Paragraph(name()).addStyle(ReportResources.chapterTitleStyle()));
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

        document.add(new Paragraph("TODO: Add other potentially other germline variants").addStyle(ReportResources.tableTitleStyle()));
    }

    private void addMVLHAnalysis(@NotNull Document document) {
        Table table = TableUtil.createReportContentTable(new float[] { 2, 2, 1, 2, 2, 1, 2, 2, 1, 2, 2, 1 },
                new Cell[] { TableUtil.createHeaderCell("Gene"), TableUtil.createHeaderCell("MVLH"), TableUtil.createHeaderCell(""),
                        TableUtil.createHeaderCell("Gene"), TableUtil.createHeaderCell("MVLH"), TableUtil.createHeaderCell(""),
                        TableUtil.createHeaderCell("Gene"), TableUtil.createHeaderCell("MVLH"), TableUtil.createHeaderCell(""),
                        TableUtil.createHeaderCell("Gene"), TableUtil.createHeaderCell("MVLH"), TableUtil.createHeaderCell("") });

        int count = 0;
        Set<String> genes = Sets.newTreeSet(Comparator.naturalOrder());
        genes.addAll(report.germlineMVLHPerGene().keySet());
        for (String gene : genes) {
            double mvlh = report.germlineMVLHPerGene().get(gene);
            if (mvlh > 0.01) {
                count++;
                table.addCell(TableUtil.createContentCell(gene));
                table.addCell(TableUtil.createContentCell(PERCENTAGE_FORMAT.format(mvlh * 100)));
                table.addCell(TableUtil.createContentCell(""));
            }
        }

        // Make sure all rows are properly filled in case table is sparse.
        if (count % 4 != 0) {
            for (int i = 0; i < 12 - 3 * (count % 4); i++) {
                table.addCell(TableUtil.createContentCell(""));
            }
        }

        document.add(TableUtil.createWrappingReportTable(table, "Genes with missed variant likelihood > 1% (" + count + ")"));
    }

    private void addGermlineCNAberrations(@NotNull Document document) {
        document.add(new Paragraph("TODO: Add Germline CN aberrations").addStyle(ReportResources.tableTitleStyle()));
    }

    private void addPEACH(@NotNull Document document) {
        document.add(new Paragraph("TODO: Add PEACH").addStyle(ReportResources.tableTitleStyle()));
    }

}
