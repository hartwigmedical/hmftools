package com.hartwig.hmftools.orange.report.tables;

import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_COHOR_FREQ;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_FUSION;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_JUNC_END;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_JUNC_START;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_SUPPORT;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_SV_TYPE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_TYPE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.cellArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSupportField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.intToFloatArray;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.isofox.RnaFusion;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Chromosomes;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

public final class RnaFusionTable
{
    private static final String COL_OTHER_FRAGS = "Non-Split Support";

    public static Table build(final String title, float width, final List<RnaFusion> fusions, final ReportResources reportResources)
    {
        if(fusions.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<Cell> cellEntries = Lists.newArrayList();

        addEntry(cells, widths, cellEntries, 2, COL_FUSION);
        addEntry(cells, widths, cellEntries, 2, COL_TYPE);
        addEntry(cells, widths, cellEntries, 1, COL_SV_TYPE);
        addEntry(cells, widths, cellEntries, 1, COL_JUNC_START);
        addEntry(cells, widths, cellEntries, 1, COL_JUNC_END);
        addEntry(cells, widths, cellEntries, 1, COL_SUPPORT);
        addEntry(cells, widths, cellEntries, 1, COL_OTHER_FRAGS);
        addEntry(cells, widths, cellEntries, 1, COL_COHOR_FREQ);

        Table table = Tables.createContent(width, intToFloatArray(widths), cellArray(cellEntries));

        for(RnaFusion fusion : sort(fusions))
        {
            table.addCell(cells.createContent(fusion.display()));

            table.addCell(cells.createContent(fusion.knownType().toString()));
            table.addCell(cells.createContent(fusion.svType().toString()));
            table.addCell(cells.createContent(fusion.junctionTypeStart()));
            table.addCell(cells.createContent(fusion.junctionTypeEnd()));

            int splitFragments = fusion.splitFragments();
            int averageDepth = min((int)round((fusion.depthStart() + fusion.depthEnd()) * 0.5), splitFragments);
            table.addCell(cells.createContent(formatSupportField(splitFragments, averageDepth)));

            table.addCell(cells.createContent(String.valueOf(fusion.realignedFrags() + fusion.discordantFrags())));
            table.addCell(cells.createContent(String.valueOf(fusion.cohortFrequency())));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    private static List<RnaFusion> sort(final List<RnaFusion> fusions)
    {
        return fusions.stream().sorted((fusion1, fusion2) ->
        {
            String locationUp1 = Chromosomes.zeroPrefixed(fusion1.chromosomeStart());
            String locationUp2 = Chromosomes.zeroPrefixed(fusion2.chromosomeStart());

            if(locationUp1.equals(locationUp2))
            {
                return Integer.compare(fusion1.positionStart(), fusion2.positionStart());
            }
            else
            {
                return locationUp1.compareTo(locationUp2);
            }
        }).collect(Collectors.toList());
    }
}
