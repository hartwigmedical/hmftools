package com.hartwig.hmftools.orange.report.chapters;

import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.PurpleQCInterpretation;
import com.hartwig.hmftools.orange.report.tables.HLAAlleleTable;
import com.hartwig.hmftools.orange.report.tables.ImmuneEscapeTable;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;

import org.jetbrains.annotations.NotNull;

public class ImmunologyChapter implements ReportChapter
{
    @NotNull
    private final OrangeRecord report;
    @NotNull
    private final ReportResources reportResources;

    public ImmunologyChapter(@NotNull final OrangeRecord report, @NotNull final ReportResources reportResources)
    {
        this.report = report;
        this.reportResources = reportResources;
    }

    @NotNull
    @Override
    public String name()
    {
        return "Immunology";
    }

    @NotNull
    @Override
    public PageSize pageSize()
    {
        return PageSize.A4;
    }

    @Override
    public void render(@NotNull final Document document)
    {
        document.add(new Paragraph(name()).addStyle(reportResources.chapterTitleStyle()));

        addHLAData(document);
        addImmuneEscapeData(document);
    }

    private void addHLAData(@NotNull Document document)
    {
        Cells cells = new Cells(reportResources);
        Table qc = new Table(UnitValue.createPercentArray(new float[] { 1, 1 }));
        qc.addCell(cells.createKey("QC Status:"));
        qc.addCell(cells.createValue(report.lilac().qc()));

        document.add(new Tables(reportResources).createWrapping(qc, "HLA QC"));

        String title = "HLA Alleles (" + report.lilac().alleles().size() + ")";
        boolean isTumorFail = PurpleQCInterpretation.isFail(report.purple().fit().qc());
        document.add(HLAAlleleTable.build(title, contentWidth(), report.lilac().alleles(), reportResources, isTumorFail));
    }

    private void addImmuneEscapeData(@NotNull Document document)
    {
        String title = "Genetic Immune Escape";
        boolean isTumorFail = PurpleQCInterpretation.isFail(report.purple().fit().qc());
        document.add(ImmuneEscapeTable.build(title, contentWidth(), report.immuneEscape(), reportResources, isTumorFail));
    }
}
