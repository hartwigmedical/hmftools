package com.hartwig.hmftools.orange.report.chapters;

import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
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
    private final OrangeRecord report;

    public ImmunologyChapter(@NotNull final OrangeRecord report) {
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
    public void render(@NotNull final Document document, @NotNull ReportResources reportResources) {
        document.add(new Paragraph(name()).addStyle(reportResources.chapterTitleStyle()));

        addHLAData(document, reportResources);
    }

    private void addHLAData(final Document document, @NotNull ReportResources reportResources) {
        Cells cells = reportResources.cells();
        Table qc = new Table(UnitValue.createPercentArray(new float[] { 1, 1 }));
        qc.addCell(cells.createKey("QC Status:"));
        qc.addCell(cells.createValue(report.lilac().qc()));

        document.add(reportResources.tables().createWrapping(qc, "HLA QC"));

        String title = "HLA Alleles (" + report.lilac().alleles().size() + ")";
        document.add(HLAAlleleTable.build(title, contentWidth(), report.lilac().alleles(), reportResources));
    }
}
