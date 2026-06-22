package com.hartwig.hmftools.orange.report.chapters;

import static com.hartwig.hmftools.common.hla.HlaCommon.MHC_CLASS_I;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.hla.LilacAllele;
import com.hartwig.hmftools.datamodel.hla.LilacRecord;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.pdfdata.ImmunologyData;
import com.hartwig.hmftools.orange.report.tables.HLAAlleleTable;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Paragraph;

import org.jetbrains.annotations.NotNull;

public class ImmunologyChapter implements ReportChapter
{
    private final ImmunologyData mData;
    private final ReportResources mReportResources;

    public ImmunologyChapter(final ImmunologyData data, final ReportResources reportResources)
    {
        mData = data;
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
    public PageSize pageSize()
    {
        return PageSize.A4;
    }

    @Override
    public void render(@NotNull final Document document)
    {
        document.add(new Paragraph(name()).addStyle(mReportResources.chapterTitleStyle()));

        if(mData.hasPurpleFail)
        {
            mReportResources.addQcFailNotice(document);
            return;
        }

        addHLAData(document);
    }

    private void addHLAData(final Document document)
    {
        LilacRecord lilacData = mData.lilac;

        if(lilacData == null)
        {
            return;
        }

        String title = "HLA Class I Alleles";

        List<LilacAllele> classIAlleles =
                lilacData.alleles().stream().filter(x -> x.geneClass().equals(MHC_CLASS_I)).collect(Collectors.toList());
        document.add(HLAAlleleTable.build(title, contentWidth(), classIAlleles, mReportResources, mData.hasRna));

        if(mData.isWholeGenome)
        {
            title = "HLA Class II Alleles";
            List<LilacAllele> classIIAlleles =
                    lilacData.alleles().stream().filter(x -> !x.geneClass().equals(MHC_CLASS_I)).collect(Collectors.toList());
            document.add(HLAAlleleTable.build(title, contentWidth(), classIIAlleles, mReportResources, mData.hasRna));
        }
    }
}
