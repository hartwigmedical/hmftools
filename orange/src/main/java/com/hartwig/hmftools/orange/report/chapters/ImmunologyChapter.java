package com.hartwig.hmftools.orange.report.chapters;

import static com.hartwig.hmftools.common.hla.HlaCommon.MHC_CLASS_I;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.hla.LilacAllele;
import com.hartwig.hmftools.datamodel.hla.LilacRecord;
import com.hartwig.hmftools.datamodel.orange.ExperimentType;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.orange.report.DocumentContext;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.algo.QcStatusInterpretation;
import com.hartwig.hmftools.orange.report.tables.HLAAlleleTable;
import com.hartwig.hmftools.orange.report.tables.ImmuneEscapeTable;

import org.apache.pdfbox.pdmodel.common.PDRectangle;
import org.jetbrains.annotations.NotNull;

public class ImmunologyChapter implements ReportChapter
{
    private final OrangeRecord mReport;
    private final ReportResources mReportResources;

    public ImmunologyChapter(final OrangeRecord report, final ReportResources reportResources)
    {
        mReport = report;
        mReportResources = reportResources;
    }

    @NotNull
    @Override
    public String name()
    {
        return "Immunology";
    }

    @NotNull
    @Override
    public PDRectangle pageSize()
    {
        return PDRectangle.A4;
    }

    @Override
    public void render(@NotNull final DocumentContext document) throws IOException
    {
        document.addParagraph(name(), mReportResources.chapterTitleStyle());

        if(QcStatusInterpretation.hasPurpleFail(mReport.purple().fit().qc()))
        {
            document.addQcFailNotice(mReportResources);
            return;
        }

        addHLAData(document);

        if(mReport.experimentType() == ExperimentType.WHOLE_GENOME)
        {
            addImmuneEscapeData(document);
        }
    }

    private void addHLAData(final DocumentContext document) throws IOException
    {
        LilacRecord lilacData = mReport.lilac();

        if(lilacData == null)
        {
            return;
        }

        String title = "HLA Class I Alleles";

        List<LilacAllele> classIAlleles =
                lilacData.alleles().stream().filter(x -> x.geneClass().equals(MHC_CLASS_I)).collect(Collectors.toList());
        document.addTable(HLAAlleleTable.build(document, title, contentWidth(), classIAlleles, mReportResources, mReport.hasRna()));

        if(mReport.experimentType() == ExperimentType.WHOLE_GENOME)
        {
            title = "HLA Class II Alleles";
            List<LilacAllele> classIIAlleles =
                    lilacData.alleles().stream().filter(x -> !x.geneClass().equals(MHC_CLASS_I)).collect(Collectors.toList());
            document.addTable(HLAAlleleTable.build(document, title, contentWidth(), classIIAlleles, mReportResources, mReport.hasRna()));
        }
    }

    private void addImmuneEscapeData(final DocumentContext document) throws IOException
    {
        String title = "Genetic Immune Escape";
        boolean isTumorFail = QcStatusInterpretation.hasPurpleFail(mReport.purple().fit().qc());
        document.addTable(ImmuneEscapeTable.build(document, title, contentWidth(), mReport.immuneEscape(), mReportResources, isTumorFail));
    }
}
