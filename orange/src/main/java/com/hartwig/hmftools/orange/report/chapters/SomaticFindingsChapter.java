package com.hartwig.hmftools.orange.report.chapters;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.sigs.SignatureAllocation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;
import com.hartwig.hmftools.datamodel.finding.DriverInterpretation;
import com.hartwig.hmftools.orange.report.PlotPathResolver;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.datamodel.finding.BreakendEntry;
import com.hartwig.hmftools.datamodel.finding.BreakendEntryFactory;
import com.hartwig.hmftools.datamodel.finding.SmallVariant;
import com.hartwig.hmftools.datamodel.finding.SmallVariantFactory;
import com.hartwig.hmftools.datamodel.finding.PurpleQCInterpretation;
import com.hartwig.hmftools.datamodel.finding.VariantDedup;
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
    @NotNull
    private final OrangeRecord report;
    @NotNull
    private final PlotPathResolver plotPathResolver;
    @NotNull
    private final ReportResources reportResources;

    public SomaticFindingsChapter(@NotNull final OrangeRecord report, @NotNull final PlotPathResolver plotPathResolver,
            @NotNull final ReportResources reportResources)
    {
        this.report = report;
        this.plotPathResolver = plotPathResolver;
        this.reportResources = reportResources;
    }

    @NotNull
    @Override
    public String name()
    {
        return "Somatic Findings";
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

        addSomaticVariants(document);
        if(!PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
        {
            addKataegisPlot(document);
        }
        addSomaticAmpDels(document);
        addFusions(document);

        if(!report.tumorOnlyMode())
        {
            addViralPresence(document);
        }

        addHomozygousDisruptions(document);
        addBreakends(document);

        if(!report.tumorOnlyMode())
        {
            addSignatureAllocations(document);
        }

        if(!PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
        {
            addStructuralDriverPlots(document);
        }
    }

    private void addSomaticVariants(@NotNull Document document)
    {
        String driverVariantsTitle = "Driver variants";
        String otherPotentiallyInterestingTitle = "Other potentially relevant variants";

        if(PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
        {
            Tables tables = new Tables(reportResources);
            document.add(tables.createNotAvailable(driverVariantsTitle, contentWidth()));
            document.add(tables.createNotAvailable(otherPotentiallyInterestingTitle, contentWidth()));
        }
        else
        {
            List<PurpleDriver> somaticDrivers = report.purple().somaticDrivers();

            List<SmallVariant> reportableVariants =
                    SmallVariantFactory.create(VariantDedup.apply(report.purple().reportableSomaticVariants()), somaticDrivers);
            String titleDrivers = driverVariantsTitle + " (" + reportableVariants.size() + ")";
            document.add(SomaticVariantTable.build(titleDrivers, contentWidth(), reportableVariants, reportResources));
        }
    }

    private void addKataegisPlot(@NotNull Document document)
    {
        document.add(new Paragraph("Kataegis plot").addStyle(reportResources.tableTitleStyle()));
        String kataegisPlot = report.plots().purpleKataegisPlot();
        if(kataegisPlot != null)
        {
            Image image = Images.build(plotPathResolver.resolve(kataegisPlot));
            image.setMaxWidth(contentWidth());
            image.setHorizontalAlignment(HorizontalAlignment.CENTER);
            document.add(image);
        }
        else
        {
            document.add(new Paragraph("No kataegis plot could be generated for this sample")
                    .addStyle(reportResources.tableContentStyle()));
        }
    }

    private void addSomaticAmpDels(@NotNull Document document)
    {
        String driverAmpsDelsTitle = "Driver amplifications and homozygous deletions";

        if(PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
        {
            Tables tables = new Tables(reportResources);
            document.add(tables.createNotAvailable(driverAmpsDelsTitle, contentWidth()));
        }
        else
        {
            String titleDrivers = driverAmpsDelsTitle + " (" + report.purple().reportableSomaticGainsDels().size() + ")";
            document.add(GainDeletionTable.build(titleDrivers,
                    contentWidth(),
                    report.purple().reportableSomaticGainsDels(),
                    report.isofox(),
                    reportResources));
        }
    }

    private void addFusions(@NotNull Document document)
    {
        String driverFusionsTitle = "Driver fusions";

        if(PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
        {
            Tables tables = new Tables(reportResources);
            document.add(tables.createNotAvailable(driverFusionsTitle, contentWidth()));
        }
        else
        {
            String titleDrivers = driverFusionsTitle + " (" + report.linx().reportableSomaticFusions().size() + ")";
            document.add(DnaFusionTable.build(titleDrivers,
                    contentWidth(),
                    report.linx().reportableSomaticFusions(),
                    report.isofox(),
                    reportResources));
        }
    }

    private void addViralPresence(@NotNull Document document)
    {
        VirusInterpreterData virusInterpreter = report.virusInterpreter();

        if(virusInterpreter != null)
        {
            String driverVirusTitle = "Driver viruses";

            if(PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
            {
                Tables tables = new Tables(reportResources);
                document.add(tables.createNotAvailable(driverVirusTitle, contentWidth()));
            }
            else
            {
                String titleDrivers = driverVirusTitle + " (" + virusInterpreter.reportableViruses().size() + ")";
                document.add(ViralPresenceTable.build(titleDrivers, contentWidth(), virusInterpreter.reportableViruses(), reportResources));
            }
        }
    }

    private void addHomozygousDisruptions(@NotNull Document document)
    {
        String homozygousDisruptionTitle = "Homozygous disruptions";

        if(PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
        {
            Tables tables = new Tables(reportResources);
            document.add(tables.createNotAvailable(homozygousDisruptionTitle, contentWidth()));
        }
        else
        {
            String title = homozygousDisruptionTitle + " (" + report.linx().somaticHomozygousDisruptions().size() + ")";
            document.add(HomozygousDisruptionTable.build(title,
                    contentWidth(),
                    report.linx().somaticHomozygousDisruptions(),
                    reportResources));
        }
    }

    private void addBreakends(@NotNull Document document)
    {
        String driverGeneDisruptionsTitle = "Driver gene disruptions";
        String nonDriverGeneDisruptionsTitle = "Other potentially interesting gene disruptions";

        if(PurpleQCInterpretation.isContaminated(report.purple().fit().qc()))
        {
            Tables tables = new Tables(reportResources);
            document.add(tables.createNotAvailable(driverGeneDisruptionsTitle, contentWidth()));
            document.add(tables.createNotAvailable(nonDriverGeneDisruptionsTitle, contentWidth()));
        }
        else
        {
            List<BreakendEntry> reportableBreakends = BreakendEntryFactory.create(report.linx().reportableSomaticBreakends(),
                    report.linx().allSomaticStructuralVariants(),
                    report.linx().somaticDrivers());

            String titleDriver = driverGeneDisruptionsTitle + " (" + reportableBreakends.size() + ")";
            document.add(BreakendTable.build(titleDriver, contentWidth(), reportableBreakends, reportResources));
        }
    }

    private void addSignatureAllocations(@NotNull Document document)
    {
        List<SignatureAllocation> sigAllocations = report.sigAllocations();

        if(sigAllocations != null)
        {
            String signatureTitle = "Signature allocations";

            if(PurpleQCInterpretation.isFail(report.purple().fit().qc()))
            {
                Tables tables = new Tables(reportResources);
                document.add(tables.createNotAvailable(signatureTitle, contentWidth()));
            }
            else
            {
                String title = signatureTitle + " (" + sigAllocations.size() + ")";
                document.add(SignatureAllocationTable.build(title, contentWidth(), sigAllocations, reportResources));
            }
        }
    }

    private void addStructuralDriverPlots(@NotNull Document document)
    {
        String title = "Structural driver plots (" + report.plots().linxDriverPlots().size() + ")";
        document.add(new Paragraph(title).addStyle(reportResources.tableTitleStyle()));
        Table table = new Table(2);
        Cells cells = new Cells(reportResources);
        for(String plot : report.plots().linxDriverPlots())
        {
            Image image = Images.build(plotPathResolver.resolve(plot));
            image.setMaxWidth(Math.round(contentWidth() / 2D) - 2);
            image.setHorizontalAlignment(HorizontalAlignment.CENTER);
            table.addCell(cells.createImage(image));
        }

        if(report.plots().linxDriverPlots().size() % 2 == 1)
        {
            table.addCell(cells.createContent(Strings.EMPTY));
        }

        document.add(table);
    }
}
