package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.layout.element.Paragraph;
import org.jetbrains.annotations.NotNull;

public class BodyText {

    @NotNull
    public static final Paragraph createBodyText(@NotNull String text) {

        return new Paragraph()
            .add(text)
            .addStyle(ReportResources.bodyTextStyle());

    }

}
