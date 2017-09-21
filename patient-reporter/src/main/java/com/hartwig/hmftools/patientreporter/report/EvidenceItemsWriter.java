package com.hartwig.hmftools.patientreporter.report;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.hyperLink;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.stl;

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.apiclients.civic.api.CivicApiWrapper;
import com.hartwig.hmftools.apiclients.civic.data.CivicEvidenceItem;
import com.hartwig.hmftools.apiclients.civic.data.CivicVariant;
import com.hartwig.hmftools.apiclients.civic.data.CivicVariantType;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.variant.Variant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.patientreporter.HmfReporterData;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.PatientReporterApplication;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.builder.component.MultiPageListBuilder;
import net.sf.dynamicreports.report.builder.component.VerticalListBuilder;
import net.sf.dynamicreports.report.builder.style.StyleBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;
import net.sf.dynamicreports.report.constant.VerticalTextAlignment;
import net.sf.dynamicreports.report.exception.DRException;

public class EvidenceItemsWriter {

    private static final Logger LOGGER = LogManager.getLogger(EvidenceItemsWriter.class);

    private static final String FONT = "Times New Roman";
    private static final Color BORKIE_COLOR = new Color(221, 235, 247);

    private static final int TEXT_HEADER_INDENT = 30;
    private static final int TEXT_DETAIL_INDENT = 40;
    private static final int LIST_INDENT = 5;
    private static final int HEADER_TO_DETAIL_VERTICAL_GAP = 8;
    private static final int DETAIL_TO_DETAIL_VERTICAL_GAP = 4;
    private static final int SECTION_VERTICAL_GAP = 25;
    private static final int PADDING = 1;

    @NotNull
    private final String reportDirectory;

    public EvidenceItemsWriter(@NotNull final String reportDirectory) {
        this.reportDirectory = reportDirectory;
    }

    //    @Override
    public void writeSequenceReport(@NotNull final PatientReport report, @NotNull final HmfReporterData reporterData)
            throws IOException, DRException {
        final JasperReportBuilder reportBuilder = generatePatientReport(report, reporterData);
        writeReport(report.sample(), reportBuilder);
    }

    private void writeReport(@NotNull final String sample, @NotNull final JasperReportBuilder report)
            throws FileNotFoundException, DRException {
        final String fileName = fileName(sample);
        if (Files.exists(new File(fileName).toPath())) {
            LOGGER.warn(" Could not write report as it already exists: " + fileName);
        } else {
            report.toPdf(new FileOutputStream(fileName));
            LOGGER.info(" Created patient report at " + fileName);
        }
    }

    @NotNull
    private String fileName(@NotNull final String sample) {
        return reportDirectory + File.separator + sample + "_evidence_report_detailed.pdf";
    }

    @VisibleForTesting
    @NotNull
    static JasperReportBuilder generatePatientReport(@NotNull final PatientReport report, @NotNull final HmfReporterData reporterData) {
        // @formatter:off

        final MultiPageListBuilder totalReport = cmp.multiPageList().add(
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("HMF Sequencing Report v" + PatientReporterApplication.VERSION + " - Civic Evidence Items").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP));
        for(final VariantReport variantReport : report.variants()){
            for(final HmfGenomeRegion region: reporterData.geneModel().hmfRegions()){
                if(region.gene().equals(variantReport.gene())){
                    final int entrezId = Integer.parseInt(region.entrezId());
                    final ComponentBuilder<?, ?> evidenceItemsPage =
                        cmp.verticalList(evidenceSection(entrezId, variantReport));
                    totalReport.add(evidenceItemsPage);
                    totalReport.newPage();
                }
            }
        }
        // @formatter:on
        CivicApiWrapper.releaseResources();
        return report().noData(totalReport);
    }

    @NotNull
    private static ComponentBuilder<?, ?> evidenceSection(final int entrezId, @NotNull final VariantReport variantReport) {
        final VerticalListBuilder section = toList("Civic evidence for variant: " + entrezId + "(" + variantReport.gene() + ")" + "\t"
                + variantReport.chromosomePosition() + "\t" + variantReport.variantField(), Lists.newArrayList());

        final Variant variant = variantReportToVariant(variantReport);
        final List<CivicVariant> civicVariants = CivicApiWrapper.getVariantsContaining(entrezId, variant).toList().blockingGet();

        final List<CivicVariant> exactMatchVariants = civicVariants.stream()
                .filter(civicVariant -> civicVariant.coordinates().equals(variant))
                .filter(civicVariant -> !civicVariant.groupedEvidenceItems().isEmpty())
                .collect(Collectors.toList());
        for (final CivicVariant civicVariant : exactMatchVariants) {
            section.add(cmp.verticalGap(HEADER_TO_DETAIL_VERTICAL_GAP), cmp.text(
                    "Exact match: " + civicVariant.name() + " " + civicVariant.coordinates() + " (" + Strings.join(
                            civicVariant.variantTypes().stream().map(CivicVariantType::name).collect(Collectors.toList()), ',') + ")")
                            .setHyperLink(hyperLink(civicVariant.summaryUrl()))
                            .setStyle(linkStyle()), cmp.verticalGap(HEADER_TO_DETAIL_VERTICAL_GAP),
                    //                    createEvidenceTable(civicVariant.evidenceItemsWithDrugs()), cmp.verticalGap(HEADER_TO_DETAIL_VERTICAL_GAP),
                    createConciseEvidenceTable(civicVariant.groupedEvidenceItems()));
        }
        final List<CivicVariant> approximateMatchVariants = civicVariants.stream()
                .filter(civicVariant -> !exactMatchVariants.contains(civicVariant))
                .filter(civicVariant -> !civicVariant.groupedEvidenceItems().isEmpty())
                .collect(Collectors.toList());
        for (final CivicVariant civicVariant : approximateMatchVariants) {
            section.add(cmp.verticalGap(HEADER_TO_DETAIL_VERTICAL_GAP), cmp.text(
                    "Approximate match: " + civicVariant.name() + " " + civicVariant.coordinates() + " (" + Strings.join(
                            civicVariant.variantTypes().stream().map(CivicVariantType::name).collect(Collectors.toList()), ',') + ")")
                            .setHyperLink(hyperLink(civicVariant.summaryUrl()))
                            .setStyle(linkStyle()), cmp.verticalGap(HEADER_TO_DETAIL_VERTICAL_GAP),
                    //                    createEvidenceTable(civicVariant.evidenceItemsWithDrugs()), cmp.verticalGap(HEADER_TO_DETAIL_VERTICAL_GAP),
                    createConciseEvidenceTable(civicVariant.groupedEvidenceItems()));
        }
        return section;
    }

    @NotNull
    private static ComponentBuilder<?, ?> createEvidenceTable(@NotNull final List<CivicEvidenceItem> evidenceItems) {
        //@formatter:off
        final int fontSize = 7;
        return cmp.subreport(
                baseTable().setColumnStyle(dataStyle().setFontSize(fontSize)).fields(EvidenceDataSource.evidenceFields())
                    .columns(
                        col.column("Level", EvidenceDataSource.LEVEL_FIELD).setFixedWidth(50),
                        col.column("Tumor Type", EvidenceDataSource.TUMOR_TYPE_FIELD).setFixedWidth(100),
                        col.column("Direction", EvidenceDataSource.DIRECTION_FIELD).setFixedWidth(100),
                        col.column("Significance", EvidenceDataSource.SIGNIFICANCE_FIELD).setFixedWidth(100),
                        col.column("Drugs", EvidenceDataSource.DRUGS_FIELD).setFixedWidth(150)))
                    .setDataSource(EvidenceDataSource.fromEvidenceItems(evidenceItems));
        // @formatter:on
    }

    @NotNull
    private static ComponentBuilder<?, ?> createConciseEvidenceTable(
            @NotNull final Map<String, Map<String, List<CivicEvidenceItem>>> groupedEvidenceItems) {
        //@formatter:off
        final int fontSize = 7;
        return cmp.subreport(
                baseTable().setColumnStyle(dataStyle().setFontSize(fontSize)).fields(EvidenceDataSource.conciseEvidenceFields())
                    .columns(
                        col.column("Significance", EvidenceDataSource.SIGNIFICANCE_FIELD).setFixedWidth(100),
                        col.column("Drugs", EvidenceDataSource.DRUGS_FIELD).setFixedWidth(400)))
                    .setDataSource(EvidenceDataSource.fromGroupedEvidenceItems(groupedEvidenceItems));
        // @formatter:on
    }

    @NotNull
    private static VerticalListBuilder toList(@NotNull final String title, @NotNull final Iterable<String> lines) {
        final VerticalListBuilder list = cmp.verticalList();
        list.add(cmp.horizontalList(cmp.horizontalGap(TEXT_HEADER_INDENT), cmp.text(title).setStyle(fontStyle().bold().setFontSize(11))),
                cmp.verticalGap(HEADER_TO_DETAIL_VERTICAL_GAP));
        boolean isFirst = true;
        for (final String line : lines) {
            if (!isFirst) {
                list.add(cmp.verticalGap(DETAIL_TO_DETAIL_VERTICAL_GAP));
            }
            list.add(cmp.horizontalList(cmp.horizontalGap(TEXT_DETAIL_INDENT), cmp.text("- ").setStyle(fontStyle()).setWidth(LIST_INDENT),
                    cmp.text(line).setStyle(fontStyle())));

            isFirst = false;
        }
        return list;
    }

    @NotNull
    private static JasperReportBuilder baseTable() {
        return report().setColumnStyle(dataStyle()).setColumnTitleStyle(tableHeaderStyle()).highlightDetailEvenRows();
    }

    @NotNull
    private static StyleBuilder sectionHeaderStyle() {
        return fontStyle().bold().setFontSize(12).setHorizontalTextAlignment(HorizontalTextAlignment.CENTER);
    }

    @NotNull
    private static StyleBuilder tableHeaderStyle() {
        return fontStyle().bold()
                .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)
                .setVerticalTextAlignment(VerticalTextAlignment.MIDDLE)
                .setFontSize(10)
                .setBorder(stl.pen1Point())
                .setBackgroundColor(BORKIE_COLOR)
                .setPadding(PADDING);
    }

    @NotNull
    private static StyleBuilder dataStyle() {
        return fontStyle().setFontSize(8)
                .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)
                .setVerticalTextAlignment(VerticalTextAlignment.MIDDLE)
                .setPadding(PADDING);
    }

    @NotNull
    private static StyleBuilder linkStyle() {
        return dataStyle().setForegroundColor(Color.BLUE);
    }

    @NotNull
    private static StyleBuilder fontStyle() {
        return stl.style().setFontName(FONT);
    }

    private static Variant variantReportToVariant(@NotNull final VariantReport variantReport) {
        return new Variant() {
            @NotNull
            @Override
            public String ref() {
                return variantReport.ref();
            }

            @NotNull
            @Override
            public String alt() {
                return variantReport.alt();
            }

            @NotNull
            @Override
            public VariantType type() {
                return null;
            }

            @NotNull
            @Override
            public String filter() {
                return null;
            }

            @NotNull
            @Override
            public String chromosome() {
                return variantReport.chromosome();
            }

            @Override
            public long position() {
                return variantReport.position();
            }
        };
    }

}
