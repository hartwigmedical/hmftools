package com.hartwig.hmftools.patientreporter.report.pages;

import static com.hartwig.hmftools.patientreporter.report.Commons.SECTION_VERTICAL_GAP;
import static com.hartwig.hmftools.patientreporter.report.Commons.dataTableStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.fontStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.formattedDate;
import static com.hartwig.hmftools.patientreporter.report.Commons.tableHeaderStyle;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;

import com.hartwig.hmftools.common.lims.LimsSampleType;
import com.hartwig.hmftools.patientreporter.QCFailReport;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailStudy;
import com.hartwig.hmftools.patientreporter.report.Commons;
import com.hartwig.hmftools.patientreporter.report.components.MainPageTopSection;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.builder.component.TextFieldBuilder;

@Value.Immutable
@Value.Style(passAnnotations = NotNull.class,
             allParameters = true)
public abstract class QCFailPage {

    @NotNull
    abstract SampleReport sampleReport();

    @NotNull
    abstract String user();

    @NotNull
    abstract QCFailReason reason();

    @NotNull
    abstract QCFailStudy study();

    @NotNull
    public static QCFailPage of(@NotNull QCFailReport report) {
        return ImmutableQCFailPage.of(report.sampleReport(), report.user(), report.reason(), report.study());
    }

    @NotNull
    public ComponentBuilder<?, ?> reportComponent() {
        return cmp.verticalList(MainPageTopSection.build(reason().title(), sampleReport()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                mainPageQCFailSection(sampleReport().sampleId()));
    }

    @NotNull
    private ComponentBuilder<?, ?> mainPageQCFailSection(@NotNull String sampleId) {
        LimsSampleType type = LimsSampleType.fromSampleId(sampleId);

        if (sampleReport().recipient() == null) {
            throw new IllegalStateException("No recipient address present for sample " + sampleReport().sampleId());
        }

        final String title;
        final String subTitle;
        final String message;

        switch (reason()) {
            case LOW_DNA_YIELD: {
                title = "Notification tumor sample on hold for sequencing";
                subTitle = "Insufficient amount of DNA";
                message = "The amount of isolated DNA was <50 ng, which is insufficient for sequencing. ";
                break;
            }
            case LOW_TUMOR_PERCENTAGE: {
                title = "Notification of inadequate tumor sample";
                subTitle = "Not enough tumor cells detected by Pathology UMC Utrecht.";
                message = "For sequencing we require a minimum of 30% tumor cells.";
                break;
            }
            case POST_ANALYSIS_FAIL: {
                title = "Notification of inadequate tumor sample";
                subTitle = "Analysis has failed post DNA isolation";
                message = "This sample could not be processed to completion successfully.";
                break;
            }
            case SHALLOW_SEQ_LOW_PURITY: {
                title = "Notification of inadequate tumor sample";
                subTitle = "Not enough tumor DNA detected by molecular T % estimate.";
                message = "For sequencing we require a minimum of 30% tumor cells.";
                break;
            }
            default: {
                title = "TITLE";
                subTitle = "SUB_TITLE";
                message = "MESSAGE";
            }
        }

        return type == LimsSampleType.CORE ? CORELayout(title, subTitle, message) : CPCTDRUPLayout(title, subTitle, message);
    }

    @NotNull
    private ComponentBuilder<?, ?> CORELayout(@NotNull String title, @NotNull String subTitle, @NotNull String message) {
        String contactDetails = sampleReport().contactNames() + " (" + sampleReport().contactEmails() + ")";
        return cmp.verticalList(cmp.text(title).setStyle(tableHeaderStyle().setFontSize(12)).setHeight(20),
                cmp.text(subTitle).setStyle(dataTableStyle().setFontSize(12)).setHeight(20),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text(message).setStyle(fontStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                notSequencedText(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("When possible, please resubmit using the same DVO with project name " + sampleReport().projectName() + ".")
                        .setStyle(fontStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("The HMF sample ID is " + sampleReport().sampleId() + " and the hospital patient ID is "
                        + sampleReport().hospitalPatientId()).setStyle(fontStyle()),
                cmp.text("The project name of sample is " + sampleReport().projectName() + " and the submission ID is "
                        + sampleReport().submission()).setStyle(fontStyle()),
                cmp.text("The internal tumor barcode is " + sampleReport().barcodeTumor() + " and the internal blood barcode is "
                        + sampleReport().barcodeReference()).setStyle(fontStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("The tumor percentage estimated by Pathology UMC Utrecht is " + sampleReport().pathologyTumorPercentage())
                        .setStyle(fontStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                shallowSeqText(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                sampleArrivalDateText(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                recipientText(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("The contact details are : " + contactDetails).setStyle(fontStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                accreditationText(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                questionsText());
    }

    @NotNull
    private ComponentBuilder<?, ?> CPCTDRUPLayout(@NotNull String title, @NotNull String subTitle, @NotNull String message) {
        return cmp.verticalList(cmp.text(title).setStyle(tableHeaderStyle().setFontSize(12)).setHeight(20),
                cmp.text(subTitle).setStyle(dataTableStyle().setFontSize(12)).setHeight(20),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text(message).setStyle(fontStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                notSequencedText(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("When possible, please resubmit using the same " + study().studyName() + "-number. "
                        + "In case additional tumor material cannot be provided, please be notified that the patient will not be "
                        + "evaluable for the " + study().studyCode() + " study.").setStyle(fontStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("The HMF sample ID is " + sampleReport().sampleId()).setStyle(fontStyle()),
                cmp.text("The internal tumor barcode is " + sampleReport().barcodeTumor() + " and the internal blood barcode is "
                        + sampleReport().barcodeReference()).setStyle(fontStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("The tumor percentage estimated by Pathology UMC Utrecht is " + sampleReport().pathologyTumorPercentage())
                        .setStyle(fontStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                shallowSeqText(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                sampleArrivalDateText(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                recipientText(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                accreditationText(),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                questionsText());
    }

    @NotNull
    private static TextFieldBuilder<String> questionsText() {
        return cmp.text("For questions, please contact us via info@hartwigmedicalfoundation.nl").setStyle(fontStyle());
    }

    @NotNull
    private TextFieldBuilder<String> shallowSeqText() {
        return cmp.text("The tumor percentage estimated by molecular tumor percentage is  " + sampleReport().purityShallowSeq())
                .setStyle(fontStyle());
    }

    @NotNull
    private TextFieldBuilder<String> sampleArrivalDateText() {
        return cmp.text(
                "The biopsies evaluated for this sample have arrived on " + formattedDate(sampleReport().tumorArrivalDate()) + " at "
                        + Commons.HARTWIG_ADDRESS).setStyle(fontStyle());
    }

    @NotNull
    private static TextFieldBuilder<String> notSequencedText() {
        return cmp.text("The received biopsies for the tumor sample for this patient were inadequate to obtain a reliable sequencing "
                + "result. Therefore whole genome sequencing cannot be performed, "
                + "unless additional fresh tumor material can be provided for a new assessment.").setStyle(fontStyle());
    }

    @NotNull
    private TextFieldBuilder<String> recipientText() {
        return cmp.text("This report is generated and verified by: " + user() + " and is addressed at " + sampleReport().recipient())
                .setStyle(fontStyle());
    }

    @NotNull
    private TextFieldBuilder<String> accreditationText() {
        return cmp.text("The results on this report are based on tests that are performed under ISO/ICE-17025:2005 accreditation.")
                .setStyle(fontStyle());
    }

}
