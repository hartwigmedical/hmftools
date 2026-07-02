package com.hartwig.hmftools.orange.report.tables;

import static java.lang.Math.round;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_COHOR_FREQ;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_FUSION;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_JUNC_END;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_JUNC_START;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_SUPPORT;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_SV_TYPE;
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
import com.hartwig.hmftools.datamodel.isofox.RnaFusion;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;

import be.quodlibet.boxable.BaseTable;

import com.hartwig.hmftools.orange.report.DocumentContext;

import java.io.IOException;

public final class RnaFusionTable
{
    public static BaseTable build(final DocumentContext docCtx, final String title, float width, final List<RnaFusion> fusions,
            final ReportResources reportResources) throws IOException
    {
        if(fusions.isEmpty())
        {
            return createEmptyTable(docCtx, title, width, reportResources);
        }

        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<String> headers = Lists.newArrayList();

        addEntry(widths, headers, 2, COL_FUSION);
        addEntry(widths, headers, 2, COL_TYPE);
        addEntry(widths, headers, 1, COL_SV_TYPE);
        addEntry(widths, headers, 1, COL_JUNC_START);
        addEntry(widths, headers, 1, COL_JUNC_END);
        addEntry(widths, headers, 1, COL_SUPPORT);
        // addEntry(widths, headers, 1, COL_OTHER_FRAGS);
        addEntry(widths, headers, 1, COL_COHOR_FREQ);

        BaseTable table = createStandardTable(docCtx, title, width, intToFloatArray(widths), stringArray(headers), reportResources);
        float[] pcts = toPercentages(intToFloatArray(widths));

        for(RnaFusion fusion : sort(fusions))
        {
            List<String> rowValues = Lists.newArrayList();
            rowValues.add(fusion.display());

            rowValues.add(fusion.knownType().toString());
            rowValues.add(fusion.svType().toString());
            rowValues.add(fusion.junctionTypeStart());
            rowValues.add(fusion.junctionTypeEnd());

            int splitFragments = fusion.splitFragments();
            int averageDepth = (int) round((fusion.depthStart() + fusion.depthEnd()) * 0.5);
            rowValues.add(formatSupportField(splitFragments, averageDepth));

            // rowValues.add(String.valueOf(fusion.realignedFrags() + fusion.discordantFrags()));
            rowValues.add(String.valueOf(fusion.cohortFrequency()));
            cells.addRow(table, pcts, rowValues);
        }

        return table;
    }

    private static List<RnaFusion> sort(final List<RnaFusion> fusions)
    {
        return fusions.stream().sorted((fusion1, fusion2) ->
        {
            String locationUp1 = zeroPrefixed(fusion1.chromosomeStart());
            String locationUp2 = zeroPrefixed(fusion2.chromosomeStart());

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
