package com.hartwig.hmftools.patientreporter.report.pages;

import static com.hartwig.hmftools.patientreporter.report.Commons.SECTION_VERTICAL_GAP;
import static com.hartwig.hmftools.patientreporter.report.Commons.dataTableStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.fontStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.formattedDate;
import static com.hartwig.hmftools.patientreporter.report.Commons.tableHeaderStyle;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;

import com.hartwig.hmftools.patientreporter.NotAnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.algo.NotAnalysableReason;
import com.hartwig.hmftools.patientreporter.algo.NotAnalysableStudy;
import com.hartwig.hmftools.patientreporter.report.Commons;
import com.hartwig.hmftools.patientreporter.report.components.MainPageTopSection;

import org.apache.logging.log4j.util.Strings;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;

@Value.Immutable
@Value.Style(passAnnotations = NotNull.class,
             allParameters = true)
public abstract class NonSequenceablePage {

    @NotNull
    abstract SampleReport sampleReport();

    @NotNull
    abstract String user();

    @NotNull
    abstract NotAnalysableReason reason();

    @NotNull
    abstract NotAnalysableStudy study();

    @NotNull
    public static NonSequenceablePage of(@NotNull final NotAnalysedPatientReport report) {
        return ImmutableNonSequenceablePage.of(report.sampleReport(), report.user(), report.reason(), report.study());
    }

    @NotNull
    public ComponentBuilder<?, ?> reportComponent() {
        return cmp.verticalList(MainPageTopSection.buildWithPathologyTumorPercentage(reason().title(), sampleReport()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                mainPageNotAnalysableSection());
    }

    @NotNull
    private ComponentBuilder<?, ?> mainPageNotAnalysableSection() {
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
            default: {
                title = "TITLE";
                subTitle = "SUB_TITLE";
                message = "MESSAGE";
            }
        }

        return cmp.verticalList(cmp.text(title).setStyle(tableHeaderStyle().setFontSize(12)).setHeight(20),
                cmp.text(subTitle).setStyle(dataTableStyle().setFontSize(12)).setHeight(20),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text(message).setStyle(fontStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("The received biopsies for the tumor sample for this patient were inadequate to obtain a reliable sequencing "
                        + "result. Therefore whole genome sequencing cannot be performed, "
                        + "unless additional fresh tumor material can be provided for a new assessment.").setStyle(fontStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("When possible, please resubmit using the same " + study().studyName() + "-number. "
                        + "In case additional tumor material cannot be provided, please be notified that the patient will not be "
                        + "evaluable for the " + study().studyCode() + " study.").setStyle(fontStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("The biopsies evaluated for this sample have arrived on " + formattedDate(sampleReport().tumorArrivalDate())
                        + " at " + Commons.HARTWIG_ADDRESS).setStyle(fontStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("This report is generated and verified by: " + user() + " and is addressed at " + sampleReport().recipient())
                        .setStyle(fontStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("The results on this report are based is generated from tests that are performed under ISO/ICE-17025:2005 accreditation.").setStyle(fontStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("For questions, please contact us via info@hartwigmedicalfoundation.nl").setStyle(fontStyle()));
    }
}
