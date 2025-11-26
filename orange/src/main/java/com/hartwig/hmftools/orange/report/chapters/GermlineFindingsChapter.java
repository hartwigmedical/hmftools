package com.hartwig.hmftools.orange.report.chapters;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.finding.GainDeletion;
import com.hartwig.hmftools.datamodel.finding.SmallVariant;
import com.hartwig.hmftools.datamodel.linx.LinxHomozygousDisruption;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.peach.PeachGenotype;
import com.hartwig.hmftools.datamodel.purple.PurpleGermlineAberration;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.finding.BreakendEntry;
import com.hartwig.hmftools.orange.report.finding.BreakendEntryFactory;
import com.hartwig.hmftools.orange.report.tables.BreakendTable;
import com.hartwig.hmftools.orange.report.tables.GainDeletionTable;
import com.hartwig.hmftools.orange.report.tables.GermlineVariantTable;
import com.hartwig.hmftools.orange.report.tables.HomozygousDisruptionTable;
import com.hartwig.hmftools.orange.report.tables.MissedVariantLikelihoodTable;
import com.hartwig.hmftools.orange.report.tables.PharmacogeneticsTable;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;

import org.jetbrains.annotations.NotNull;

public class GermlineFindingsChapter implements ReportChapter
{
    @NotNull
    private final OrangeRecord report;
    @NotNull
    private final ReportResources reportResources;

    public GermlineFindingsChapter(@NotNull final OrangeRecord report, @NotNull final ReportResources reportResources)
    {
        this.report = report;
        this.reportResources = reportResources;
    }

    @NotNull
    @Override
    public String name()
    {
        return "Germline Findings";
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

        if(report.refSample() != null)
        {
            addGermlineVariants(document);
            addGermlineDeletions(document);
            addGermlineHomozygousDisruptions(document);
            addGermlineBreakends(document);
            addMVLHAnalysis(document);
            addGermlineCNAberrations(document);
            addPharmacogenetics(document);
        }
        else
        {
            document.add(new Paragraph(ReportResources.NOT_AVAILABLE).addStyle(reportResources.tableContentStyle()));
        }
    }

    private void addGermlineVariants(@NotNull Document document)
    {
        List<SmallVariant> reportableEntries = report.findings().driverGermlineSmallVariants(ReportedStatus.REPORTED);
        if(reportableEntries != null)
        {
            String titleDrivers = "Driver variants (" + reportableEntries.size() + ")";
            document.add(GermlineVariantTable.build(titleDrivers, contentWidth(), reportableEntries, reportResources));
        }
    }

    private void addGermlineDeletions(@NotNull Document document)
    {
        List<GainDeletion> reportableGermlineFullDels = report.purple().driverGermlineDeletions();
        if(reportableGermlineFullDels != null)
        {
            String title = "Potentially pathogenic germline deletions (" + reportableGermlineFullDels.size() + ")";
            document.add(GainDeletionTable.build(title, contentWidth(), reportableGermlineFullDels, report.isofox(), reportResources));
        }
    }

    private void addGermlineHomozygousDisruptions(@NotNull Document document)
    {
        List<LinxHomozygousDisruption> germlineHomozygousDisruptions = report.linx().germlineHomozygousDisruptions();
        if(germlineHomozygousDisruptions != null)
        {
            String title = "Potentially pathogenic germline homozygous disruptions (" + germlineHomozygousDisruptions.size() + ")";
            document.add(HomozygousDisruptionTable.build(title, contentWidth(), germlineHomozygousDisruptions, reportResources));
        }
    }

    private void addGermlineBreakends(@NotNull Document document)
    {
        List<LinxSvAnnotation> allGermlineStructuralVariants = report.linx().allGermlineStructuralVariants();
        List<LinxBreakend> reportableGermlineBreakends = report.linx().driverGermlineBreakends();

        if(allGermlineStructuralVariants != null && reportableGermlineBreakends != null)
        {
            // TODO: Load Linx germline drivers properly
            List<BreakendEntry> reportableBreakends =
                    BreakendEntryFactory.create(reportableGermlineBreakends, allGermlineStructuralVariants, List.of());

            String title = "Potentially pathogenic germline gene disruptions (" + reportableBreakends.size() + ")";
            document.add(BreakendTable.build(title, contentWidth(), reportableBreakends, reportResources));
        }
    }

    private void addMVLHAnalysis(@NotNull Document document)
    {
        Map<String, Double> germlineMVLHPerGene = report.germlineMVLHPerGene();
        if(germlineMVLHPerGene != null)
        {
            Map<String, Double> significantGermlineMVLHPerGene = germlineMVLHPerGene.entrySet()
                    .stream()
                    .filter(e -> e.getValue() > 0.01)
                    .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));

            String title = "Genes with missed variant likelihood > 1% (" + significantGermlineMVLHPerGene.size() + ")";
            document.add(MissedVariantLikelihoodTable.build(title, contentWidth(), significantGermlineMVLHPerGene, reportResources));
        }
    }

    private void addGermlineCNAberrations(@NotNull Document document)
    {
        Set<PurpleGermlineAberration> germlineAberrations = report.purple().fit().qc().germlineAberrations();
        if(!germlineAberrations.isEmpty())
        {
            int count = 0;
            StringJoiner germlineAberrationJoiner = new StringJoiner(", ");
            for(PurpleGermlineAberration germlineAberration : germlineAberrations)
            {
                if(germlineAberration != PurpleGermlineAberration.NONE)
                {
                    count++;
                }
                germlineAberrationJoiner.add(germlineAberration.toString());
            }
            Table table = new Table(UnitValue.createPercentArray(new float[] { 1 })).setWidth(contentWidth());
            table.addCell(new Cells(reportResources).createContent(germlineAberrationJoiner.toString()));
            document.add(new Tables(reportResources).createWrapping(table, "Germline CN aberrations (" + count + ")"));
        }
    }

    private void addPharmacogenetics(@NotNull Document document)
    {
        Set<PeachGenotype> peach = report.peach();
        if(peach != null)
        {
            String titlePharmacogenetics = "Pharmacogenetics (" + peach.size() + ")";
            document.add(PharmacogeneticsTable.build(titlePharmacogenetics, contentWidth(), peach, reportResources));
        }
    }
}
