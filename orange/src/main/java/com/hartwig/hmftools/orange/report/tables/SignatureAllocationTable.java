package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;
import static com.hartwig.hmftools.orange.report.ReportResources.formatPercentage;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.cellArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.intToFloatArray;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.sigs.SignatureAllocation;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;

public final class SignatureAllocationTable
{
    static final String MISALLOC_SIGNATURE = "MISALLOC";

    public static Table build(
            final String title, float width, final List<SignatureAllocation> signatureAllocations, final ReportResources reportResources)
    {
        if(signatureAllocations.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<Cell> cellEntries = Lists.newArrayList();

        addEntry(cells, widths, cellEntries, 1, "Signature");
        addEntry(cells, widths, cellEntries, 2, "Etiology");
        addEntry(cells, widths, cellEntries, 1, "Allocation");
        addEntry(cells, widths, cellEntries, 1, "Percent");
        addEntry(cells, widths, cellEntries, 4, Strings.EMPTY); // to space things out

        /*
        Table table = Tables.createContent(width,
                new float[] { 10, 23, 10, 10, 10 },
                new Cell[] { cells.createHeader("Signature"), cells.createHeader("Etiology"), cells.createHeader("Allocation"),
                        cells.createHeader("Percent"),
                        cells.createHeader(Strings.EMPTY) });
        */

        Table table = Tables.createContent(width, intToFloatArray(widths), cellArray(cellEntries));

        for(SignatureAllocation signatureAllocation : sort(signatureAllocations))
        {
            table.addCell(cells.createContent(signatureAllocation.signature()));
            table.addCell(cells.createContent(signatureAllocation.etiology()));
            table.addCell(cells.createContent(format("%.0f", signatureAllocation.allocation())));
            table.addCell(cells.createContent(formatPercentage(signatureAllocation.percent())));
            table.addCell(cells.createContent(Strings.EMPTY));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    @VisibleForTesting
    static List<SignatureAllocation> sort(final List<SignatureAllocation> signatureAllocations)
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

            return Double.compare(allocation2.allocation(), allocation1.allocation());
        }).collect(Collectors.toList());
    }
}
