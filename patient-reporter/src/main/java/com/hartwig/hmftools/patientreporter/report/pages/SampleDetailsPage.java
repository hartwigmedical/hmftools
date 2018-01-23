package com.hartwig.hmftools.patientreporter.report.pages;

import static com.hartwig.hmftools.patientreporter.report.Commons.SECTION_VERTICAL_GAP;
import static com.hartwig.hmftools.patientreporter.report.Commons.formattedDate;
import static com.hartwig.hmftools.patientreporter.report.Commons.sectionHeaderStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.toList;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;

import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.report.Commons;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;

@Value.Immutable
@Value.Style(passAnnotations = NotNull.class,
             allParameters = true)
public abstract class SampleDetailsPage {

    @NotNull
    abstract SampleReport sampleReport();

    @NotNull
    abstract String user();

    @NotNull
    abstract Optional<String> comments();

    public static SampleDetailsPage of(@NotNull final PatientReport report) {
        return ImmutableSampleDetailsPage.of(report.sampleReport(), report.user(), report.comments());
    }

    public ComponentBuilder<?, ?> reportComponent() {
        return cmp.verticalList(cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text(Commons.TITLE + " - Sample Details & Disclaimer").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP), sampleDetailsSection(), cmp.verticalGap(SECTION_VERTICAL_GAP), disclaimerSection());
    }

    @NotNull
    private static ComponentBuilder<?, ?> disclaimerSection() {
        //@formatter:off
        final List<String> lines = Lists.newArrayList(
                "The data on which this report is based is generated from tests that are performed under ISO/ICE-17025:2005 accreditation.",
                "This analysis done for this report has passed all internal quality controls.",
                "For feedback or complaints please contact qualitysystem@hartwigmedicalfoundation.nl.");
        //@formatter:on
        return toList("Disclaimer", lines);
    }

    @NotNull
    private ComponentBuilder<?, ?> sampleDetailsSection() {
        if (sampleReport().recipient() == null) {
            throw new IllegalStateException("No recipient address present for sample " + sampleReport().sampleId());
        }

        //@formatter:off
        final List<String> lines = Lists.newArrayList(
                "The samples have been sequenced at Hartwig Medical Foundation, Science Park 408, 1098XH Amsterdam",
                "The samples have been analyzed by Next Generation Sequencing",
                "This experiment is performed on the tumor sample which arrived on " +
                        formattedDate(sampleReport().tumorArrivalDate()),
                "The pathology tumor percentage for this sample is " + sampleReport().tumorPercentageString(),
                "This experiment is performed on the blood sample which arrived on " +
                        formattedDate(sampleReport().bloodArrivalDate()),
                "This experiment is performed according to lab procedures: " + sampleReport().labProcedures(),
                "This report is generated and verified by: " + user(),
                "This report is addressed at: " + sampleReport().recipient());
        //@formatter:on
        comments().ifPresent(comments -> lines.add("Comments: " + comments));

        return toList("Sample details", lines);
    }
}
