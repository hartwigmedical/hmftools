package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.layout.element.Paragraph;
import org.jetbrains.annotations.NotNull;

public class BodyText {

    public static final Paragraph getParagraph(@NotNull String text) {
        Paragraph p = new Paragraph();
        p.add(text);
        p.addStyle(ReportResources.bodyTextStyle());
        return p;
    }

}
