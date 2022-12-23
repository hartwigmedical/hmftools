package com.hartwig.hmftools.orange.report.chapters;

import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.tables.HLAAlleleTable;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;

import org.jetbrains.annotations.NotNull;

public class ImmunologyChapter implements ReportChapter {

    @NotNull
    private final OrangeReport report;

    public ImmunologyChapter(@NotNull final OrangeReport report) {
        this.report = report;
    }

    @NotNull
    @Override
    public String name() {
        return "Immunology";
    }

    @NotNull
    @Override
    public PageSize pageSize() {
        return PageSize.A4;
    }

    @Override
    public void render(@NotNull final Document document) {
        document.add(new Paragraph(name()).addStyle(ReportResources.chapterTitleStyle()));

        addHLAData(document);
    }

    private void addHLAData(final Document document) {
        Table qc = new Table(UnitValue.createPercentArray(new float[] { 1, 1 }));
        qc.addCell(Cells.createKey("QC Status:"));
        qc.addCell(Cells.createValue(report.lilac().qc()));

        document.add(Tables.createWrapping(qc, "HLA QC"));

        String title = "HLA Alleles (" + report.lilac().alleles().size() + ")";
        document.add(HLAAlleleTable.build(title, contentWidth(), report.lilac().alleles()));
    }
}
