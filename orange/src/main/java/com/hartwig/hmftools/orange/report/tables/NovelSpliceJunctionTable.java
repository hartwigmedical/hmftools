package com.hartwig.hmftools.orange.report.tables;

import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_COHOR_FREQ;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_GENE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_JUNCTIONS;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_JUNC_END;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_JUNC_START;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_SUPPORT;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_TYPE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.stringArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createStandardTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createEmptyTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.toPercentages;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSupportField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.intToFloatArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.zeroPrefixed;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.isofox.AltSpliceJunctionType;
import com.hartwig.hmftools.datamodel.isofox.NovelSpliceJunction;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;

import be.quodlibet.boxable.BaseTable;

import com.hartwig.hmftools.orange.report.DocumentContext;

import java.io.IOException;

public final class NovelSpliceJunctionTable
{
    public static BaseTable build(final DocumentContext docCtx,
            final String title, float width, final List<NovelSpliceJunction> junctions, final ReportResources reportResources)
            throws IOException
    {
        if(junctions.isEmpty())
        {
            return createEmptyTable(docCtx, title, width, reportResources);
        }

        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<String> headers = Lists.newArrayList();

        addEntry(widths, headers, 1, COL_GENE);
        addEntry(widths, headers, 3, COL_JUNCTIONS);
        addEntry(widths, headers, 2, COL_TYPE);
        addEntry(widths, headers, 1, COL_JUNC_START);
        addEntry(widths, headers, 1, COL_JUNC_END);
        addEntry(widths, headers, 1, COL_SUPPORT);
        addEntry(widths, headers, 1, COL_COHOR_FREQ);

        BaseTable table = createStandardTable(docCtx, title, width, intToFloatArray(widths), stringArray(headers), reportResources);
        float[] pcts = toPercentages(intToFloatArray(widths));

        for(NovelSpliceJunction junction : sort(junctions))
        {
            List<String> rowValues = Lists.newArrayList();
            rowValues.add(junction.gene());
            rowValues.add(junctionsDisplay(junction));
            rowValues.add(junction.type().toString());
            rowValues.add(String.valueOf(junction.regionStart()));
            rowValues.add(String.valueOf(junction.regionEnd()));

            int fragments = junction.fragmentCount();
            int averageDepth = (int) round((junction.depthStart() + junction.depthEnd()) * 0.5);
            rowValues.add(formatSupportField(fragments, averageDepth));

            rowValues.add(String.valueOf(junction.cohortFrequency()));
            cells.addRow(table, pcts, rowValues);
        }

        return table;
    }

    private static String junctionsDisplay(final NovelSpliceJunction junction)
    {
        boolean dupType = junction.type() == AltSpliceJunctionType.CIRCULAR;

        int positionStart = dupType ? junction.junctionEnd() : junction.junctionStart();
        int positionEnd = dupType ? junction.junctionStart() : junction.junctionEnd();
        int exonStart = dupType ? junction.exonEnd() : junction.exonStart();
        int exonEnd = dupType ? junction.exonStart() : junction.exonEnd();

        return format("Exon %d (%s:%d) - Exon %d (%s:%d)",
                exonStart, junction.chromosome(), positionStart, exonEnd, junction.chromosome(), positionEnd);
    }

    private static List<NovelSpliceJunction> sort(final List<NovelSpliceJunction> junctions)
    {
        return junctions.stream().sorted((junction1, junction2) ->
        {
            String locationUp1 = zeroPrefixed(junction1.chromosome());
            String locationUp2 = zeroPrefixed(junction2.chromosome());

            if(locationUp1.equals(locationUp2))
            {
                return Integer.compare(junction1.junctionStart(), junction2.junctionStart());
            }
            else
            {
                return locationUp1.compareTo(locationUp2);
            }
        }).collect(Collectors.toList());
    }
}
