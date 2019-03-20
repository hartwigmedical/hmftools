package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.patientreporter.cfreport.components.BodyText;
import com.hartwig.hmftools.patientreporter.cfreport.components.LineDivider;
import com.hartwig.hmftools.patientreporter.cfreport.components.SectionTitle;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import org.jetbrains.annotations.NotNull;

public class SummaryChapter extends ReportChapter {

    @Override
    public final String getName() {
        return "Summary";
    }

    @Override
    public final ChapterType getChapterType() {
        return ChapterType.SummaryChapter;
    }

    @Override
    protected final void renderChapterContent(@NotNull Document report) {
        renderTumorLocationAndType(report);
        renderSummaryText(report);
        renderTreatmentIndications(report);
        renderTumorCharacteristics(report);
        renderGenomicAlterations(report);
    }

    private final void renderTumorLocationAndType(@NotNull Document report) {

        Div div = new Div();
        div.setKeepTogether(true);
        div.setWidth(getContentWidth());

    }

    private final void renderSummaryText(@NotNull Document report) {

        Div div = initializeSummarySectionDiv("Summary", getContentWidth());

        // Add content to div
        String content = "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Vivamus eget porta turpis. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nullam interdum sodales ullamcorper. Nulla vestibulum ipsum quis ipsum congue, quis commodo velit condimentum. Suspendisse eget nulla egestas, fermentum urna ut, bibendum ipsum. Nulla varius, dui elementum faucibus ultricies, nisi velit dignissim arcu, nec feugiat dui magna eu felis. Maecenas at odio pharetra, sodales velit vitae, gravida mauris. Pellentesque id ultrices diam. Integer non ex ut neque auctor pellentesque. Ut et nibh faucibus, pretium erat efficitur, vehicula lorem. Ut in fermentum velit, et aliquet mi. Praesent et commodo lorem. Nulla eu sem nec purus tempus auctor. Vivamus vitae varius dui, a fermentum nulla. Sed scelerisque sollicitudin blandit. Ut vitae augue pretium, ultrices purus id, porttitor arcu. Cras et tortor non diam interdum venenatis. In at dolor et odio placerat feugiat id quis ante.";
        Paragraph p = BodyText.getParagraph(content);
        p.setWidth(getContentWidth());
        div.add(p);

        report.add(div);

    }

    private final void renderTreatmentIndications(@NotNull Document report) {

        Div div = initializeSummarySectionDiv("Treatment indications", getContentWidth());

        // @TODO add content

        report.add(div);

    }

    private final void renderTumorCharacteristics(@NotNull Document report) {

        Div div = initializeSummarySectionDiv("Tumor characteristics summary", getContentWidth());

        // @TODO add content

        report.add(div);

    }

    private final void renderGenomicAlterations(@NotNull Document report) {

        Div div = initializeSummarySectionDiv("Genomic alterations summary", getContentWidth());

        // @TODO add content

        report.add(div);

    }

    /**
     *
     * @param sectionTitle
     * @return
     */
    private static final Div initializeSummarySectionDiv(@NotNull String sectionTitle, float width) {

        Div div = new Div();
        div.setKeepTogether(true);
        div.setWidth(width);

        // Add divider and section title
        div.add(LineDivider.getLineDivider(width));
        div.add(SectionTitle.getSectionTitle(sectionTitle));

        return div;

    }

}
