package com.hartwig.hmftools.orange.report.chapters;

import static com.hartwig.hmftools.common.hla.HlaCommon.MHC_CLASS_I;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.hla.LilacAllele;
import com.hartwig.hmftools.datamodel.hla.LilacRecord;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.datamodel.purple.PurpleQCInterpretation;
import com.hartwig.hmftools.orange.report.tables.HLAAlleleTable;
import com.hartwig.hmftools.orange.report.tables.ImmuneEscapeTable;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;

public class ImmunologyChapter implements ReportChapter
{
    private final OrangeRecord mReport;
    private final ReportResources mReportResources;

    public ImmunologyChapter(final OrangeRecord report, final ReportResources reportResources)
    {
        mReport = report;
        mReportResources = reportResources;
    }

    @Override
    public String name()
    {
        return "Immunology";
    }

    @Override
    public PageSize pageSize()
    {
        return PageSize.A4;
    }

    @Override
    public void render(final Document document)
    {
        document.add(new Paragraph(name()).addStyle(mReportResources.chapterTitleStyle()));

        addHLAData(document);
        addImmuneEscapeData(document);
    }

    private void addHLAData(final Document document)
    {
        LilacRecord lilacData = mReport.lilac();

        if(lilacData == null)
            return;

        // Cells cells = new Cells(mReportResources);
        // Table qc = new Table(UnitValue.createPercentArray(new float[] { 1, 1 }));

        String title = "HLA Class I Alleles";

        List<LilacAllele> classIAlleles = lilacData.alleles().stream().filter(x -> x.geneClass().equals(MHC_CLASS_I)).collect(Collectors.toList());
        document.add(HLAAlleleTable.build(title, contentWidth(), classIAlleles, mReportResources, mReport.hasRna()));

        title = "HLA Class II Alleles";
        List<LilacAllele> classIIAlleles = lilacData.alleles().stream().filter(x -> !x.geneClass().equals(MHC_CLASS_I)).collect(Collectors.toList());
        document.add(HLAAlleleTable.build(title, contentWidth(), classIIAlleles, mReportResources, mReport.hasRna()));
    }

    private void addImmuneEscapeData(final Document document)
    {
        String title = "Genetic Immune Escape";
        boolean isTumorFail = PurpleQCInterpretation.isFail(mReport.purple().fit().qc());
        document.add(ImmuneEscapeTable.build(title, contentWidth(), mReport.immuneEscape(), mReportResources, isTumorFail));
    }
}
