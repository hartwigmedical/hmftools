package com.hartwig.hmftools.orange.report.util;

import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class DocumentUtil {

    private DocumentUtil() {
    }

    public static void addCheckedTable(@NotNull Document document, @NotNull String title, @Nullable Table table) {
        if (table == null) {
            document.add(TableUtil.createNoneReportTable(title));
        } else {
            document.add(table);
        }
    }
}
