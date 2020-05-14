package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.common.lims.LimsStudy;
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
        PdfCanvas canvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());
        Rectangle pageSize = page.getPageSize();
        renderBackgroundRect(fullHeight, canvas, pageSize);
        BaseMarker.renderMarkerGrid(4, (fullHeight ? 20 : 2), CONTENT_X_START, 35, 820, -ROW_SPACING, .05f, .15f, canvas);

        int sideTextIndex = -1;
        Canvas cv = new Canvas(canvas, page.getDocument(), page.getPageSize());

        cv.add(createSidePanelDiv(++sideTextIndex, "HMF sample id", sampleReport.tumorSampleId()));
        cv.add(createSidePanelDiv(++sideTextIndex, "Report date", ReportResources.REPORT_DATE));

        LimsStudy study = LimsStudy.fromSampleId(sampleReport.tumorSampleId());

        if (fullHeight && fullContent) {
            String contactNames =
                    study == LimsStudy.CORE || study == LimsStudy.WIDE ? sampleReport.hospitalData().requesterName() : Strings.EMPTY;
            if (!contactNames.isEmpty()) {
                cv.add(createSidePanelDiv(++sideTextIndex, "Requested by", contactNames));
            }

            final String contactEmails =
                    study == LimsStudy.CORE || study == LimsStudy.WIDE ? sampleReport.hospitalData().requesterEmail() : Strings.EMPTY;
            if (!contactEmails.isEmpty()) {
                cv.add(createSidePanelDiv(++sideTextIndex, "Email", contactEmails));
            }

            String hospitalName = sampleReport.hospitalData().hospitalName();
            if (!hospitalName.isEmpty()) {
                cv.add(createSidePanelDiv(++sideTextIndex, "Hospital", hospitalName));
            }

            String hospitalPatientId = study == LimsStudy.CORE ? sampleReport.hospitalPatientId() : Strings.EMPTY;
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
        float Y_START = 802;
        float VALUE_TEXT_Y_OFFSET = 18;
        float MAX_WIDTH = 120;

        Div div = new Div();
        div.setKeepTogether(true);

        float yPos = Y_START - index * ROW_SPACING;
        div.add(new Paragraph(label.toUpperCase()).addStyle(ReportResources.sidePanelLabelStyle())
                .setFixedPosition(CONTENT_X_START, yPos, MAX_WIDTH));

        float valueFontSize = ReportResources.maxPointSizeForWidth(ReportResources.fontBold(), 11, 6, value, MAX_WIDTH);
        yPos -= VALUE_TEXT_Y_OFFSET;
        div.add(new Paragraph(value).addStyle(ReportResources.sidePanelValueStyle().setFontSize(valueFontSize))
                .setHeight(15)
                .setFixedPosition(CONTENT_X_START, yPos, MAX_WIDTH)
                .setFixedLeading(valueFontSize));

        return div;
    }
}
