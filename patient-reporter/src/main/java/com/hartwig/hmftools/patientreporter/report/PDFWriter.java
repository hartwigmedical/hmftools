package com.hartwig.hmftools.patientreporter.report;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.stl;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.NotAnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.SequencedReportData;
import com.hartwig.hmftools.patientreporter.report.pages.ImmutableCircosPage;
import com.hartwig.hmftools.patientreporter.report.pages.ImmutableEvidencePage;
import com.hartwig.hmftools.patientreporter.report.pages.ImmutableExplanationPage;
import com.hartwig.hmftools.patientreporter.report.pages.ImmutableFindingsPage;
import com.hartwig.hmftools.patientreporter.report.pages.NonSequenceablePage;
import com.hartwig.hmftools.patientreporter.report.pages.SampleDetailsPage;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.dynamicreports.report.exception.DRException;

public class PDFWriter {

    private static final Logger LOGGER = LogManager.getLogger(PDFWriter.class);

    @NotNull
    private final String reportDirectory;

    public PDFWriter(@NotNull final String reportDirectory) {
        this.reportDirectory = reportDirectory;
    }

    public void writeSequenceReport(@NotNull final AnalysedPatientReport report, @NotNull final SequencedReportData reporterData)
            throws IOException, DRException {
        final JasperReportBuilder reportBuilder = generatePatientReport(report, reporterData);
        writeReport(fileName(report.sampleReport().sampleId()), reportBuilder);
    }

    public void writeNonSequenceableReport(@NotNull final NotAnalysedPatientReport report) throws IOException, DRException {
        final JasperReportBuilder reportBuilder = generateNotAnalysableReport(report);
        writeReport(fileName(report.sampleReport().sampleId()), reportBuilder);
    }

    private static void writeReport(@NotNull final String fileName, @NotNull final JasperReportBuilder report)
            throws FileNotFoundException, DRException {
        if (Files.exists(new File(fileName).toPath())) {
            LOGGER.warn(" Could not write " + fileName + " as it already exists.");
        } else {
            report.toPdf(new FileOutputStream(fileName));
            LOGGER.info(" Created patient report at " + fileName);
        }
    }

    @NotNull
    private String fileName(@NotNull final String sample) {
        return reportDirectory + File.separator + sample + "_hmf_report.pdf";
    }

    @VisibleForTesting
    @NotNull
    static JasperReportBuilder generateNotAnalysableReport(@NotNull final NotAnalysedPatientReport report) {
        // MIVO: hack to get page footers working; the footer band and noData bands are exclusive, see additional comment below for details
        final DRDataSource singleItemDataSource = new DRDataSource("item");
        singleItemDataSource.add(new Object());

        return report().pageFooter(cmp.pageXslashY())
                .lastPageFooter(cmp.verticalList(logoRVAfooter(report.logoRVA())))
                .lastPageFooter(cmp.verticalList(signatureFooter(report.signaturePath()),
                        cmp.pageXslashY(),
                        cmp.text("End of report.").setStyle(stl.style().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER))))
                .addDetail(NonSequenceablePage.of(report).reportComponent())
                .setDataSource(singleItemDataSource);
    }

    @VisibleForTesting
    @NotNull
    static JasperReportBuilder generatePatientReport(@NotNull final AnalysedPatientReport report,
            @NotNull final SequencedReportData reporterData) {
        final ComponentBuilder<?, ?> totalReport = cmp.multiPageList()
                .add(ImmutableEvidencePage.of(report).reportComponent())
                .newPage()
                .add(ImmutableFindingsPage.of(report).reportComponent())
                .newPage()
                .add(ImmutableCircosPage.of(report.circosPath()).reportComponent())
                .newPage()
                .add(ImmutableExplanationPage.builder().build().reportComponent())
                .newPage()
                .add(SampleDetailsPage.of(report).reportComponent());

        // MIVO: hack to get page footers working; the footer band and noData bands are exclusive:
        //  - footerBand, detailBand, etc are shown when data source is not empty
        //  - noData band is shown when data source is empty; intended to be used when there is no data to show in the report
        //  (e.g. would be appropriate to be used for notAnalysableReport)
        //
        // more info: http://www.dynamicreports.org/examples/bandreport

        final DRDataSource singleItemDataSource = new DRDataSource("item");
        singleItemDataSource.add(new Object());

        return report().pageFooter(cmp.pageXslashY())
                .lastPageFooter(cmp.verticalList(logoRVAfooter(report.logoRVA())))
                .lastPageFooter(cmp.verticalList(signatureFooter(report.signaturePath()),
                        cmp.pageXslashY(),
                        cmp.text("End of report.").setStyle(stl.style().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER))))
                .addDetail(totalReport)
                .setDataSource(singleItemDataSource);
    }

    @NotNull
    private static ComponentBuilder<?, ?> logoRVAfooter(@NotNull final String logoPath) {
        return cmp.horizontalList(cmp.horizontalGap(370), cmp.xyList().add(40, 5, cmp.image(logoPath)), cmp.verticalGap(120));
    }

    @NotNull
    private static ComponentBuilder<?, ?> signatureFooter(@NotNull final String signaturePath) {
        return cmp.horizontalList(cmp.horizontalGap(370),
                cmp.xyList()
                        .add(40, 5, cmp.image(signaturePath))
                        .add(0, 0, cmp.text("Edwin Cuppen,"))
                        .add(0, 15, cmp.text("Director Hartwig Medical Foundation").setWidth(190)),
                cmp.horizontalGap(10));
    }
}
