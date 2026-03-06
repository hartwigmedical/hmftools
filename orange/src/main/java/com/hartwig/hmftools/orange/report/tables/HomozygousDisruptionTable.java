package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_GENE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_LOCATION;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_TYPE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.cellArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.intToFloatArray;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.linx.LinxHomozygousDisruption;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Chromosomes;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;

public final class HomozygousDisruptionTable
{
    public static Table build(
            final String title, float width, final List<LinxHomozygousDisruption> homozygousDisruptions,
            final ReportResources reportResources)
    {
        if(homozygousDisruptions.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<Cell> cellEntries = Lists.newArrayList();

        addEntry(cells, widths, cellEntries, 1, COL_LOCATION);
        addEntry(cells, widths, cellEntries, 1, COL_GENE);
        addEntry(cells, widths, cellEntries, 1, COL_TYPE);
        addEntry(cells, widths, cellEntries, 1, COL_DRIVER);
        addEntry(cells, widths, cellEntries, 3, Strings.EMPTY); // to space things out

        Table table = Tables.createContent(width, intToFloatArray(widths), cellArray(cellEntries));

        for(LinxHomozygousDisruption homozygousDisruption : sort(homozygousDisruptions))
        {
            table.addCell(cells.createContent(homozygousDisruption.chromosome() + homozygousDisruption.chromosomeBand()));
            table.addCell(cells.createContent(gene(homozygousDisruption)));
            table.addCell(cells.createContent(homozygousDisruption.type()));
            table.addCell(cells.createContent(homozygousDisruption.driverInterpretation().toString()));
            table.addCell(cells.createContent(Strings.EMPTY));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    private static List<LinxHomozygousDisruption> sort(final List<LinxHomozygousDisruption> homozygousDisruptions)
    {
        return homozygousDisruptions.stream().sorted((disruption1, disruption2) ->
        {
            String location1 = Chromosomes.zeroPrefixed(disruption1.chromosome() + disruption1.chromosomeBand());
            String location2 = Chromosomes.zeroPrefixed(disruption2.chromosome() + disruption2.chromosomeBand());

            if(location1.equals(location2))
            {
                return disruption1.gene().compareTo(disruption2.gene());
            }
            else
            {
                return location1.compareTo(location2);
            }
        }).collect(Collectors.toList());
    }

    private static String gene(final LinxHomozygousDisruption homozygousDisruption)
    {
        String addon = Strings.EMPTY;
        if(!homozygousDisruption.isCanonical())
        {
            addon = " (alt)";
        }
        return homozygousDisruption.gene() + addon;
    }
}
