package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;
import static com.hartwig.hmftools.orange.report.ReportResources.formatPercentage;
import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.datamodel.sigs.SignatureAllocation;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class SignatureAllocationTable
{
    static final String MISALLOC_SIGNATURE = "MISALLOC";

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<SignatureAllocation> signatureAllocations,
            @NotNull ReportResources reportResources)
    {
        if(signatureAllocations.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);
        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 3 },
                new Cell[] { cells.createHeader("Signature"), cells.createHeader("Allocation"), cells.createHeader("Percent"),
                        cells.createHeader(Strings.EMPTY) });

        for(SignatureAllocation signatureAllocation : sort(signatureAllocations))
        {
            table.addCell(cells.createContent(signatureAllocation.signature()));
            table.addCell(cells.createContent(formatSingleDigitDecimal(signatureAllocation.allocation())));
            table.addCell(cells.createContent(formatPercentage(signatureAllocation.percent())));
            table.addCell(cells.createContent(Strings.EMPTY));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    @NotNull
    @VisibleForTesting
    static List<SignatureAllocation> sort(@NotNull List<SignatureAllocation> signatureAllocations)
    {
        return signatureAllocations.stream().sorted((allocation1, allocation2) ->
        {
            String sig1 = allocation1.signature().toLowerCase();
            String sig2 = allocation2.signature().toLowerCase();

            if(sig1.equalsIgnoreCase(MISALLOC_SIGNATURE))
            {
                return 1;
            }
            else if(sig2.equalsIgnoreCase(MISALLOC_SIGNATURE))
            {
                return -1;
            }

            if(!sig1.startsWith("sig") || !sig2.startsWith("sig"))
            {
                LOGGER.warn("Signatures do not start with 'sig': {} & {}", sig1, sig2);
                return sig1.compareTo(sig2);
            }

            int sigNumber1 = Integer.parseInt(sig1.substring(3));
            int sigNumber2 = Integer.parseInt(sig2.substring(3));

            return Integer.compare(sigNumber1, sigNumber2);
        }).collect(Collectors.toList());
    }
}
