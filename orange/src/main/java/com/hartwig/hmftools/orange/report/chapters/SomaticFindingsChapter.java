package com.hartwig.hmftools.orange.report.chapters;

import java.util.List;

import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleQCInterpretation;
import com.hartwig.hmftools.datamodel.sigs.SignatureAllocation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.orange.report.PlotPathResolver;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.datamodel.BreakendEntry;
import com.hartwig.hmftools.orange.report.datamodel.BreakendEntryFactory;
import com.hartwig.hmftools.orange.report.datamodel.VariantEntry;
import com.hartwig.hmftools.orange.report.datamodel.VariantEntryFactory;
import com.hartwig.hmftools.orange.report.interpretation.VariantDedup;
import com.hartwig.hmftools.orange.report.tables.BreakendTable;
import com.hartwig.hmftools.orange.report.tables.DnaFusionTable;
import com.hartwig.hmftools.orange.report.tables.GainDeletionTable;
import com.hartwig.hmftools.orange.report.tables.HomozygousDisruptionTable;
import com.hartwig.hmftools.orange.report.tables.SignatureAllocationTable;
import com.hartwig.hmftools.orange.report.tables.SomaticVariantTable;
import com.hartwig.hmftools.orange.report.tables.ViralPresenceTable;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Images;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.HorizontalAlignment;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class SomaticFindingsChapter implements ReportChapter
{
    private final OrangeRecord mReport;
    private final PlotPathResolver mPlotPathResolver;
    private final ReportResources mReportResources;

    public SomaticFindingsChapter(final OrangeRecord report, final PlotPathResolver plotPathResolver,
            final ReportResources reportResources)
    {
        mReport = report;
        mPlotPathResolver = plotPathResolver;
        mReportResources = reportResources;
    }

    @Override
    public String name()
    {
        return "Somatic Findings";
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

        addSomaticVariants(document);

        addSomaticAmpDels(document);
        addFusions(document);

        if(!mReport.tumorOnlyMode())
        {
            addViralPresence(document);
        }

        addHomozygousDisruptions(document);
        addBreakends(document);

        if(!mReport.tumorOnlyMode())
        {
            addSignatureAllocations(document);
        }

        if(!PurpleQCInterpretation.isContaminated(mReport.purple().fit().qc()))
        {
            addStructuralDriverPlots(document);
        }

        if(!PurpleQCInterpretation.isContaminated(mReport.purple().fit().qc()))
        {
            addKataegisPlot(document);
        }
    }

    private void addSomaticVariants(final Document document)
    {
        String driverVariantsTitle = "Driver variants";

        if(PurpleQCInterpretation.isContaminated(mReport.purple().fit().qc()))
        {
            Tables tables = new Tables(mReportResources);
            document.add(tables.createNotAvailable(driverVariantsTitle, contentWidth()));
        }
        else
        {
            List<PurpleDriver> somaticDrivers = mReport.purple().somaticDrivers();

            // List<VariantEntry> reportableVariants = VariantEntryFactory.create(VariantDedup.apply(mReport.purple().somaticVariants()), somaticDrivers);

            List<VariantEntry> reportableVariants = VariantEntryFactory.create(mReport.purple().somaticVariants(), somaticDrivers);

            String titleDrivers = driverVariantsTitle + " (" + reportableVariants.size() + ")";
            document.add(SomaticVariantTable.build(titleDrivers, contentWidth(), reportableVariants, mReportResources));
        }
    }

    private void addKataegisPlot(final Document document)
    {
        document.add(new Paragraph("Kataegis plot").addStyle(mReportResources.tableTitleStyle()));
        String kataegisPlot = mReport.plots().purpleKataegisPlot();
        if(kataegisPlot != null)
        {
            Image image = Images.build(mPlotPathResolver.resolve(kataegisPlot));
            image.setMaxWidth(contentWidth());
            image.setHorizontalAlignment(HorizontalAlignment.CENTER);
            document.add(image);
        }
        else
        {
            document.add(new Paragraph("No kataegis plot could be generated for this sample")
                    .addStyle(mReportResources.tableContentStyle()));
        }
    }

    private void addSomaticAmpDels(final Document document)
    {
        String driverAmpsDelsTitle = "Driver amplifications and homozygous deletions";

        if(PurpleQCInterpretation.isContaminated(mReport.purple().fit().qc()))
        {
            Tables tables = new Tables(mReportResources);
            document.add(tables.createNotAvailable(driverAmpsDelsTitle, contentWidth()));
        }
        else
        {
            String titleDrivers = driverAmpsDelsTitle + " (" + mReport.purple().somaticGainsDels().size() + ")";
            document.add(GainDeletionTable.build(titleDrivers,
                    contentWidth(),
                    mReport.purple().somaticGainsDels(),
                    mReport.isofox(),
                    mReportResources));
        }
    }

    private void addFusions(final Document document)
    {
        String driverFusionsTitle = "Driver fusions";

        if(PurpleQCInterpretation.isContaminated(mReport.purple().fit().qc()))
        {
            Tables tables = new Tables(mReportResources);
            document.add(tables.createNotAvailable(driverFusionsTitle, contentWidth()));
        }
        else
        {
            String titleDrivers = driverFusionsTitle + " (" + mReport.linx().fusions().size() + ")";
            document.add(DnaFusionTable.build(titleDrivers,
                    contentWidth(),
                    mReport.linx().fusions(),
                    mReport.isofox(),
                    mReportResources));
        }
    }

    private void addViralPresence(final Document document)
    {
        VirusInterpreterData virusInterpreter = mReport.virusInterpreter();

        if(virusInterpreter != null)
        {
            String driverVirusTitle = "Driver viruses";

            if(PurpleQCInterpretation.isContaminated(mReport.purple().fit().qc()))
            {
                Tables tables = new Tables(mReportResources);
                document.add(tables.createNotAvailable(driverVirusTitle, contentWidth()));
            }
            else
            {
                String titleDrivers = driverVirusTitle + " (" + virusInterpreter.reportableViruses().size() + ")";
                document.add(ViralPresenceTable.build(titleDrivers, contentWidth(), virusInterpreter.reportableViruses(), mReportResources));
            }
        }
    }

    private void addHomozygousDisruptions(final Document document)
    {
        String homozygousDisruptionTitle = "Homozygous disruptions";

        if(PurpleQCInterpretation.isContaminated(mReport.purple().fit().qc()))
        {
            Tables tables = new Tables(mReportResources);
            document.add(tables.createNotAvailable(homozygousDisruptionTitle, contentWidth()));
        }
        else
        {
            String title = homozygousDisruptionTitle + " (" + mReport.linx().somaticHomozygousDisruptions().size() + ")";
            document.add(HomozygousDisruptionTable.build(title,
                    contentWidth(),
                    mReport.linx().somaticHomozygousDisruptions(),
                    mReportResources));
        }
    }

    private void addBreakends(final Document document)
    {
        String driverGeneDisruptionsTitle = "Driver gene disruptions";
        String nonDriverGeneDisruptionsTitle = "Other potentially interesting gene disruptions";

        if(PurpleQCInterpretation.isContaminated(mReport.purple().fit().qc()))
        {
            Tables tables = new Tables(mReportResources);
            document.add(tables.createNotAvailable(driverGeneDisruptionsTitle, contentWidth()));
            document.add(tables.createNotAvailable(nonDriverGeneDisruptionsTitle, contentWidth()));
        }
        else
        {
            List<BreakendEntry> reportableBreakends = BreakendEntryFactory.create(
                    mReport.linx().somaticBreakends(),
                    mReport.linx().somaticStructuralVariants(),
                    mReport.linx().somaticDrivers());

            String titleDriver = driverGeneDisruptionsTitle + " (" + reportableBreakends.size() + ")";
            document.add(BreakendTable.build(titleDriver, contentWidth(), reportableBreakends, mReportResources));
        }
    }

    private void addSignatureAllocations(final Document document)
    {
        List<SignatureAllocation> sigAllocations = mReport.sigAllocations();

        if(sigAllocations != null)
        {
            String signatureTitle = "Signature allocations";

            if(PurpleQCInterpretation.isFail(mReport.purple().fit().qc()))
            {
                Tables tables = new Tables(mReportResources);
                document.add(tables.createNotAvailable(signatureTitle, contentWidth()));
            }
            else
            {
                String title = signatureTitle + " (" + sigAllocations.size() + ")";
                document.add(SignatureAllocationTable.build(title, contentWidth(), sigAllocations, mReportResources));
            }
        }
    }

    private void addStructuralDriverPlots(final Document document)
    {
        String title = "Structural driver plots (" + mReport.plots().linxDriverPlots().size() + ")";
        document.add(new Paragraph(title).addStyle(mReportResources.tableTitleStyle()));
        Table table = new Table(2);
        Cells cells = new Cells(mReportResources);
        for(String plot : mReport.plots().linxDriverPlots())
        {
            Image image = Images.build(mPlotPathResolver.resolve(plot));
            image.setMaxWidth(Math.round(contentWidth() / 2D) - 2);
            image.setHorizontalAlignment(HorizontalAlignment.CENTER);
            table.addCell(cells.createImage(image));
        }

        if(mReport.plots().linxDriverPlots().size() % 2 == 1)
        {
            table.addCell(cells.createContent(Strings.EMPTY));
        }

        document.add(table);
    }
}
