package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.common.lims.LimsSampleType;
import com.hartwig.hmftools.patientreporter.PatientReporterApplication;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.kernel.geom.Rectangle;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.layout.Canvas;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class SidePanel {

    private static final float ROW_SPACING = 42;
    private static final float CONTENT_X_START = 455;
    private static final float RECTANGLE_WIDTH = 170;
    private static final float RECTANGLE_HEIGHT_SHORT = 110;

    public static void renderSidePanel(PdfPage page, @NotNull final SampleReport sampleReport, boolean fullHeight, boolean fullContent) {
        final PdfCanvas canvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());
        final Rectangle pageSize = page.getPageSize();
        renderBackgroundRect(fullHeight, canvas, pageSize);
        BaseMarker.renderMarkerGrid(4, (fullHeight ? 20 : 2), CONTENT_X_START, 35, 820, -ROW_SPACING, .05f, .15f, canvas);

        int sideTextIndex = -1;
        Canvas cv = new Canvas(canvas, page.getDocument(), page.getPageSize());

        cv.add(createSidePanelDiv(++sideTextIndex, "HMF sample id", sampleReport.tumorSampleId()));
        cv.add(createSidePanelDiv(++sideTextIndex, "Report date", ReportResources.REPORT_DATE));

        LimsSampleType type = LimsSampleType.fromSampleId(sampleReport.tumorSampleId());

        if (fullHeight && fullContent) {
            final String contactNames = type == LimsSampleType.CORE
                    ? sampleReport.requesterName()
                    : type == LimsSampleType.WIDE ? sampleReport.studyRequesterName() : Strings.EMPTY;
            if (!contactNames.isEmpty()) {
                cv.add(createSidePanelDiv(++sideTextIndex, "Requested by", contactNames));
            }

            final String contactEmails = type == LimsSampleType.CORE
                    ? sampleReport.requesterEmail()
                    : type == LimsSampleType.WIDE ? sampleReport.studyRequesterEmail() : Strings.EMPTY;
            if (!contactEmails.isEmpty()) {
                cv.add(createSidePanelDiv(++sideTextIndex, "Email", contactEmails));
            }

            final String hospitalName = sampleReport.hospitalName();
            if (!hospitalName.isEmpty()) {
                cv.add(createSidePanelDiv(++sideTextIndex, "Hospital", hospitalName));
            }

            final String hospitalPatientId = type == LimsSampleType.CORE ? sampleReport.hospitalPatientId() : Strings.EMPTY;
            if (!hospitalPatientId.isEmpty()) {
                cv.add(createSidePanelDiv(++sideTextIndex, "Hospital patient id", hospitalPatientId));
            }
        }

        if (page.getDocument().getNumberOfPages() == 1) {
            cv.add(new Paragraph(
                    "v" + (PatientReporterApplication.VERSION != null ? PatientReporterApplication.VERSION : "X.X")).setFixedPosition(
                    pageSize.getWidth() - RECTANGLE_WIDTH + 4,
                    40,
                    30).setRotationAngle(Math.PI / 2).setFontColor(ReportResources.PALETTE_LIGHT_GREY).setFontSize(6));
        }

        canvas.release();
    }

    private static void renderBackgroundRect(boolean fullHeight, @NotNull PdfCanvas canvas, @NotNull Rectangle pageSize) {
        canvas.rectangle(pageSize.getWidth(),
                pageSize.getHeight(),
                -RECTANGLE_WIDTH,
                fullHeight ? -pageSize.getHeight() : -RECTANGLE_HEIGHT_SHORT);
        canvas.setFillColor(ReportResources.PALETTE_BLUE);
        canvas.fill();
    }

    @NotNull
    private static Div createSidePanelDiv(int index, @NotNull String label, @NotNull String value) {
        final float Y_START = 802;
        final float VALUE_TEXT_Y_OFFSET = 18;
        final float MAX_WIDTH = 120;

        Div div = new Div();
        div.setKeepTogether(true);

        float yPos = Y_START - index * ROW_SPACING;
        div.add(new Paragraph(label.toUpperCase()).addStyle(ReportResources.sidePanelLabelStyle())
                .setFixedPosition(CONTENT_X_START, yPos, MAX_WIDTH));

        final float valueFontSize = ReportResources.maxPointSizeForWidth(ReportResources.fontBold(), 11, 6, value, MAX_WIDTH);
        yPos -= VALUE_TEXT_Y_OFFSET;
        div.add(new Paragraph(value).addStyle(ReportResources.sidePanelValueStyle().setFontSize(valueFontSize))
                .setHeight(15)
                .setFixedPosition(CONTENT_X_START, yPos, MAX_WIDTH)
                .setFixedLeading(valueFontSize));

        return div;
    }
}
