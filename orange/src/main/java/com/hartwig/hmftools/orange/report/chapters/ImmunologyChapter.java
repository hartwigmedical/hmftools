package com.hartwig.hmftools.orange.report.chapters;

import static com.hartwig.hmftools.common.hla.HlaCommon.MHC_CLASS_I;
import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.hla.LilacAllele;
import com.hartwig.hmftools.datamodel.hla.LilacRecord;
import com.hartwig.hmftools.orange.report.OrangeFonts;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.pdfdata.ImmunologyData;
import com.hartwig.hmftools.orange.report.tables.HLAAlleleTable;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.component.VerticalListBuilder;
import net.sf.dynamicreports.report.constant.PageOrientation;
import net.sf.dynamicreports.report.constant.PageType;

public class ImmunologyChapter implements ReportChapter
{
    private final ImmunologyData mData;

    public ImmunologyChapter(final ImmunologyData data, final Object unused)
    {
        mData = data;
    }

    @Override
    public String name()
    {
        return "Immunology";
    }

    @Override
    public boolean isLandscape()
    {
        return false;
    }

    @Override
    public JasperReportBuilder buildReport()
    {
        JasperReportBuilder report = report().setPageFormat(PageType.A4, PageOrientation.PORTRAIT);
        VerticalListBuilder content = cmp.verticalList();
        content.add(cmp.text(name()).setStyle(OrangeFonts.CHAPTER_TITLE_STYLE));

        if(mData.hasPurpleFail)
        {
            content.add(cmp.text(ReportResources.NOT_AVAILABLE).setStyle(OrangeFonts.TABLE_CONTENT_STYLE));
            return report.summary(content);
        }

        LilacRecord lilac = mData.lilac;
        if(lilac == null)
        {
            return report.summary(content);
        }

        List<LilacAllele> classIAlleles = lilac.alleles().stream()
                .filter(x -> x.geneClass().equals(MHC_CLASS_I))
                .collect(Collectors.toList());
        content.add(cmp.subreport(HLAAlleleTable.build("HLA Class I Alleles", classIAlleles, null, mData.hasRna)));

        if(mData.isWholeGenome)
        {
            List<LilacAllele> classIIAlleles = lilac.alleles().stream()
                    .filter(x -> !x.geneClass().equals(MHC_CLASS_I))
                    .collect(Collectors.toList());
            content.add(cmp.subreport(HLAAlleleTable.build("HLA Class II Alleles", classIIAlleles, null, mData.hasRna)));
        }

        return report.summary(content);
    }
}
