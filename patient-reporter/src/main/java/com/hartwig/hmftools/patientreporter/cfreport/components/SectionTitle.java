package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.layout.element.Paragraph;
import org.jetbrains.annotations.NotNull;

public class SectionTitle {

    @NotNull
    public static Paragraph createSectionTitle(@NotNull String title) {

        return new Paragraph(title)
                .addStyle(ReportResources.sectionTitleStyle());

    }

}
