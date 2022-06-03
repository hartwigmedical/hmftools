package com.hartwig.hmftools.orange.report.chapters;

import java.text.DecimalFormat;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.algo.selection.GermlineVariantSelector;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.tables.GermlineDeletionTable;
import com.hartwig.hmftools.orange.report.tables.GermlineDisruptionTable;
import com.hartwig.hmftools.orange.report.tables.GermlineVariantTable;
import com.hartwig.hmftools.orange.report.tables.PharmacogeneticsTable;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class GermlineFindingsChapter implements ReportChapter {

    private static final DecimalFormat PERCENTAGE_FORMAT = ReportResources.decimalFormat("#.0'%'");

    @NotNull
    private final OrangeReport report;
    private final boolean reportGermline;

    public GermlineFindingsChapter(@NotNull final OrangeReport report, final boolean reportGermline) {
        this.report = report;
        this.reportGermline = reportGermline;
    }

    @NotNull
    @Override
    public String name() {
        return "Germline Findings";
    }

    @NotNull
    @Override
    public PageSize pageSize() {
        return PageSize.A4;
    }

    @Override
    public void render(@NotNull final Document document) {
        document.add(new Paragraph(name()).addStyle(ReportResources.chapterTitleStyle()));

        if (reportGermline) {
            addGermlineVariants(document);
            addGermlineDeletions(document);
            addGermlineDisruptions(document);
            addMVLHAnalysis(document);
            addGermlineCNAberrations(document);
            addPharmacogenetics(document);
        } else {
            document.add(new Paragraph(ReportResources.NOT_AVAILABLE).addStyle(ReportResources.tableContentStyle()));
        }
    }

    private void addGermlineVariants(@NotNull Document document) {
        String titleDrivers = "Driver variants (" + report.purple().reportableGermlineVariants().size() + ")";
        document.add(GermlineVariantTable.build(titleDrivers, contentWidth(), report.purple().reportableGermlineVariants()));

        List<ReportableVariant> nonDriverVariants = GermlineVariantSelector.selectNonDrivers(report.purple().unreportedGermlineVariants());
        String titleNonDrivers = "Other potentially relevant variants (" + nonDriverVariants.size() + ")";
        document.add(GermlineVariantTable.build(titleNonDrivers, contentWidth(), nonDriverVariants));
    }

    private void addGermlineDeletions(@NotNull Document document) {
        String title = "Potentially pathogenic germline deletions (" + report.purple().reportableGermlineDeletions().size() + ")";
        document.add(GermlineDeletionTable.build(title, contentWidth(), report.purple().reportableGermlineDeletions()));
    }

    private void addGermlineDisruptions(@NotNull Document document) {
        String title = "Potentially pathogenic germline disruptions (" + report.linx().reportableGermlineDisruptions().size() + ")";
        document.add(GermlineDisruptionTable.build(title, contentWidth(), report.linx().reportableGermlineDisruptions()));
    }

    private void addMVLHAnalysis(@NotNull Document document) {
        Table table = Tables.createContent(contentWidth(),
                new float[] { 2, 2, 1, 2, 2, 1, 2, 2, 1, 2, 2, 1 },
                new Cell[] { Cells.createHeader("Gene"), Cells.createHeader("MVLH"), Cells.createHeader(Strings.EMPTY),
                        Cells.createHeader("Gene"), Cells.createHeader("MVLH"), Cells.createHeader(Strings.EMPTY),
                        Cells.createHeader("Gene"), Cells.createHeader("MVLH"), Cells.createHeader(Strings.EMPTY),
                        Cells.createHeader("Gene"), Cells.createHeader("MVLH"), Cells.createHeader(Strings.EMPTY) });

        int count = 0;
        Set<String> genes = Sets.newTreeSet(Comparator.naturalOrder());
        genes.addAll(report.germlineMVLHPerGene().keySet());
        for (String gene : genes) {
            double mvlh = report.germlineMVLHPerGene().get(gene);
            if (mvlh > 0.01) {
                count++;
                table.addCell(Cells.createContent(gene));
                table.addCell(Cells.createContent(PERCENTAGE_FORMAT.format(mvlh * 100)));
                table.addCell(Cells.createContent(Strings.EMPTY));
            }
        }

        // Make sure all rows are properly filled in case table is sparse.
        if (count % 4 != 0) {
            for (int i = 0; i < 12 - 3 * (count % 4); i++) {
                table.addCell(Cells.createContent(Strings.EMPTY));
            }
        }

        String title = "Genes with missed variant likelihood > 1% (" + count + ")";
        if (count == 0) {
            document.add(Tables.createEmpty(title, contentWidth()));
        } else {
            document.add(Tables.createWrapping(table, title));
        }
    }

    private void addGermlineCNAberrations(@NotNull Document document) {
        int count = 0;
        StringJoiner germlineAberrations = new StringJoiner(", ");
        for (GermlineAberration aberration : report.purple().qc().germlineAberrations()) {
            if (aberration != GermlineAberration.NONE) {
                count++;
            }
            germlineAberrations.add(aberration.toString());
        }
        Table table = new Table(UnitValue.createPercentArray(new float[] { 1 })).setWidth(contentWidth());
        table.addCell(Cells.createContent(germlineAberrations.toString()));
        document.add(Tables.createWrapping(table, "Germline CN aberrations (" + count + ")"));
    }

    private void addPharmacogenetics(@NotNull Document document) {
        String titlePharmacogenetics = "Pharmacogenetics (" + report.peach().size() + ")";
        document.add(PharmacogeneticsTable.build(titlePharmacogenetics, contentWidth(), report.peach()));
    }
}
