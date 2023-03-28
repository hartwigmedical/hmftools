package com.hartwig.hmftools.orange.report.tables;

import java.text.DecimalFormat;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.datamodel.sigs.SignatureAllocation;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class SignatureAllocationTable {

    private static final Logger LOGGER = LogManager.getLogger(SignatureAllocationTable.class);

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#.#");
    private static final DecimalFormat PERCENTAGE = ReportResources.decimalFormat("#'%'");

    static final String MISALLOC_SIGNATURE = "MISALLOC";

    private SignatureAllocationTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<SignatureAllocation> signatureAllocations,
            @NotNull ReportResources reportResources) {
        if (signatureAllocations.isEmpty()) {
            return reportResources.tables().createEmpty(title, width);
        }

        Cells cells = reportResources.cells();
        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 3 },
                new Cell[] { cells.createHeader("Signature"), cells.createHeader("Allocation"), cells.createHeader("Percent"),
                        cells.createHeader(Strings.EMPTY) });

        for (SignatureAllocation signatureAllocation : sort(signatureAllocations)) {
            table.addCell(cells.createContent(signatureAllocation.signature()));
            table.addCell(cells.createContent(SINGLE_DIGIT.format(signatureAllocation.allocation())));
            table.addCell(cells.createContent(PERCENTAGE.format(signatureAllocation.percent() * 100)));
            table.addCell(cells.createContent(Strings.EMPTY));
        }

        return reportResources.tables().createWrapping(table, title);
    }

    @NotNull
    @VisibleForTesting
    static List<SignatureAllocation> sort(@NotNull List<SignatureAllocation> signatureAllocations) {
        return signatureAllocations.stream().sorted((allocation1, allocation2) -> {
            String sig1 = allocation1.signature().toLowerCase();
            String sig2 = allocation2.signature().toLowerCase();

            if (sig1.equalsIgnoreCase(MISALLOC_SIGNATURE)) {
                return 1;
            } else if (sig2.equalsIgnoreCase(MISALLOC_SIGNATURE)) {
                return -1;
            }

            if (!sig1.startsWith("sig") || !sig2.startsWith("sig")) {
                LOGGER.warn("Signatures do not start with 'sig': {} & {}", sig1, sig2);
                return sig1.compareTo(sig2);
            }

            int sigNumber1 = Integer.parseInt(sig1.substring(3));
            int sigNumber2 = Integer.parseInt(sig2.substring(3));

            return Integer.compare(sigNumber1, sigNumber2);
        }).collect(Collectors.toList());
    }
}
