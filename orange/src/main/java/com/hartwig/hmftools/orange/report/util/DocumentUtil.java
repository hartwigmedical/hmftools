package com.hartwig.hmftools.orange.report.util;

import com.hartwig.hmftools.orange.report.ReportResources;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class DocumentUtil {

    private DocumentUtil() {
    }

    public static void addCheckedTable(@NotNull Document document, @Nullable Table table, @NotNull String messageIfNull) {
        if (table == null) {
            document.add(new Paragraph(" * " + messageIfNull).addStyle(ReportResources.tableContentStyle()));
        } else {
            document.add(table);
        }
    }
}
