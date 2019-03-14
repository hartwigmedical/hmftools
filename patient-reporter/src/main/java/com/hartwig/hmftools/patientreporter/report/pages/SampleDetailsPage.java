package com.hartwig.hmftools.patientreporter.report.pages;

import static com.hartwig.hmftools.patientreporter.report.Commons.SECTION_VERTICAL_GAP;
import static com.hartwig.hmftools.patientreporter.report.Commons.formattedDate;
import static com.hartwig.hmftools.patientreporter.report.Commons.sectionHeaderStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.toList;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;

import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.lims.LimsSampleType;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.report.Commons;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;

@Value.Immutable
@Value.Style(passAnnotations = NotNull.class,
             allParameters = true)
public abstract class SampleDetailsPage {

    private static final Logger LOGGER = LogManager.getLogger(SampleDetailsPage.class);

    @NotNull
    abstract SampleReport sampleReport();

    @NotNull
    abstract String user();

    @NotNull
    abstract Optional<String> comments();

    @NotNull
    public static SampleDetailsPage of(@NotNull final PatientReport report) {
        return ImmutableSampleDetailsPage.of(report.sampleReport(), report.user(), report.comments());
    }

    @NotNull
    public ComponentBuilder<?, ?> reportComponent(@NotNull String sampleId) {
        return cmp.verticalList(cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text(Commons.TITLE_SEQUENCE + " - Sample Details & Disclaimer").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                sampleDetailsSection(sampleId),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                disclaimerSection());
    }

    @NotNull
    private static ComponentBuilder<?, ?> disclaimerSection() {
        final List<String> lines = Lists.newArrayList(
                "The data on which this report is based is generated from tests that are performed under ISO/ICE-17025:2005 accreditation.",
                "The analysis done for this report has passed all internal quality controls.",
                "For feedback or complaints please contact qualitysystem@hartwigmedicalfoundation.nl.",
                "For general questions, please contact us at info@hartwigmedicalfoundation.nl");
        return toList("Disclaimer", lines);
    }

    @NotNull
    private ComponentBuilder<?, ?> sampleDetailsSection(@NotNull String sampleId) {
        String addressee = sampleReport().addressee();
        if (addressee == null) {
            LOGGER.warn("No addressee present for sample " + sampleReport().sampleId());
        }

        addressee = addressee != null ? addressee : "N/A";

        final List<String> lines = Lists.newArrayList("The samples have been sequenced at " + Commons.HARTWIG_ADDRESS,
                "The samples have been analyzed by Next Generation Sequencing",
                "The HMF sample ID is: " + sampleReport().sampleId(),
                "The pathology tumor percentage for this sample is " + sampleReport().pathologyTumorPercentage(),
                "This experiment is performed on the tumor sample which arrived on " + formattedDate(sampleReport().tumorArrivalDate())
                        + " with internal tumor barcode " + sampleReport().barcodeTumor(),
                "This experiment is performed on the blood sample which arrived on " + formattedDate(sampleReport().bloodArrivalDate())
                        + " with internal blood barcode " + sampleReport().barcodeReference(),
                "This experiment is performed according to lab procedures: " + sampleReport().labProcedures(),
                "This report is generated and verified by: " + user(),
                "This report is addressed at: " + addressee);

        LimsSampleType type = LimsSampleType.fromSampleId(sampleId);
        if (type == LimsSampleType.CORE) {
            lines.add("The hospital patient ID is: " + sampleReport().hospitalPatientId());
            lines.add("The project name of sample is: " + sampleReport().projectName() + " and the submission ID is "
                    + sampleReport().submissionId());
            lines.add("The requester is: " + sampleReport().requesterName() + " (" + sampleReport().requesterEmail() + ")");
        }

        comments().ifPresent(comments -> lines.add("Comments: " + comments));

        return toList("Sample details", lines);
    }
}
