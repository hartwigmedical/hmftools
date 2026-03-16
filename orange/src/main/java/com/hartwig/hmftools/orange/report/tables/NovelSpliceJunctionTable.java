package com.hartwig.hmftools.orange.report.tables;

import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_COHOR_FREQ;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_GENE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_JUNC_END;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_JUNC_START;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_LOCATION;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_SUPPORT;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_TYPE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.cellArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSupportField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.intToFloatArray;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.isofox.NovelSpliceJunction;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Chromosomes;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

public final class NovelSpliceJunctionTable
{
    public static Table build(
            final String title, float width, final List<NovelSpliceJunction> junctions, final ReportResources reportResources)
    {
        if(junctions.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<Cell> cellEntries = Lists.newArrayList();

        addEntry(cells, widths, cellEntries, 1, COL_GENE);
        addEntry(cells, widths, cellEntries, 2, COL_LOCATION);
        addEntry(cells, widths, cellEntries, 2, COL_TYPE);
        addEntry(cells, widths, cellEntries, 1, COL_JUNC_START);
        addEntry(cells, widths, cellEntries, 1, COL_JUNC_END);
        addEntry(cells, widths, cellEntries, 1, COL_SUPPORT);
        addEntry(cells, widths, cellEntries, 1, COL_COHOR_FREQ);

        Table table = Tables.createContent(width, intToFloatArray(widths), cellArray(cellEntries));

        for(NovelSpliceJunction junction : sort(junctions))
        {
            table.addCell(cells.createContent(junction.gene()));
            table.addCell(cells.createContent(format("%s:%d-%d", junction.chromosome(), junction.junctionStart(), junction.junctionEnd())));
            table.addCell(cells.createContent(junction.type().toString()));
            table.addCell(cells.createContent(String.valueOf(junction.regionStart())));
            table.addCell(cells.createContent(String.valueOf(junction.regionEnd())));

            int fragments = junction.fragmentCount();
            int averageDepth = min((int)round((junction.depthStart() + junction.depthEnd()) * 0.5), fragments);
            table.addCell(cells.createContent(formatSupportField(fragments, averageDepth)));

            table.addCell(cells.createContent(String.valueOf(junction.cohortFrequency())));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    private static List<NovelSpliceJunction> sort(final List<NovelSpliceJunction> junctions)
    {
        return junctions.stream().sorted((junction1, junction2) ->
        {
            String locationUp1 = Chromosomes.zeroPrefixed(junction1.chromosome());
            String locationUp2 = Chromosomes.zeroPrefixed(junction2.chromosome());

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
