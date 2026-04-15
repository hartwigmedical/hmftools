package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatPercentage;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.stringArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createStandardTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createEmptyTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.toPercentages;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.intToFloatArray;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.sigs.SignatureAllocation;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;

import be.quodlibet.boxable.BaseTable;

import com.hartwig.hmftools.orange.report.DocumentContext;

import java.io.IOException;

import org.apache.logging.log4j.util.Strings;

public final class SignatureAllocationTable
{
    static final String MISALLOC_SIGNATURE = "MISALLOC";

    public static BaseTable build(final DocumentContext docCtx,
            final String title, float width, final List<SignatureAllocation> signatureAllocations, final ReportResources reportResources)
            throws IOException
    {
        if(signatureAllocations.isEmpty())
        {
            return createEmptyTable(docCtx, title, width, reportResources);
        }

        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<String> headers = Lists.newArrayList();

        addEntry(widths, headers, 1, "Signature");
        addEntry(widths, headers, 2, "Etiology");
        addEntry(widths, headers, 1, "Allocation");
        addEntry(widths, headers, 1, "Percent");
        addEntry(widths, headers, 4, Strings.EMPTY); // to space things out

        BaseTable table = createStandardTable(docCtx, title, width, intToFloatArray(widths), stringArray(headers), reportResources);
        float[] pcts = toPercentages(intToFloatArray(widths));

        for(SignatureAllocation signatureAllocation : sort(signatureAllocations))
        {
            List<String> rowValues = Lists.newArrayList();
            if(signatureAllocation.percent() < 0.01)
            {
                continue;
            }

            rowValues.add(signatureAllocation.signature());
            rowValues.add(signatureAllocation.etiology());
            rowValues.add(format("%.0f", signatureAllocation.allocation()));
            rowValues.add(formatPercentage(signatureAllocation.percent()));
            rowValues.add(Strings.EMPTY);
            cells.addRow(table, pcts, rowValues);
        }

        return table;
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
