package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.formatFoldChangeField;
import static com.hartwig.hmftools.orange.report.ReportResources.formatPercentileField;
import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.ReportResources.formatTpmField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_GENE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_LOCATION;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_RANGE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_REL_CN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_TPM;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_TYPE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.cellArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.intToFloatArray;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Chromosomes;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;

public final class GainDeletionTable
{
    public static Table build(
            final String title, float width, final List<PurpleGainDeletion> gainsDels, final ReportResources reportResources,
            boolean hasRna)
    {
        if(gainsDels.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<Cell> cellEntries = Lists.newArrayList();

        addEntry(cells, widths, cellEntries, 1, COL_LOCATION);
        addEntry(cells, widths, cellEntries, 1, COL_GENE);
        addEntry(cells, widths, cellEntries, 1, COL_TYPE);
        addEntry(cells, widths, cellEntries, 1, COL_RANGE); // if known
        addEntry(cells, widths, cellEntries, 1, "Min CN");
        addEntry(cells, widths, cellEntries, 1, "Max CN");
        addEntry(cells, widths, cellEntries, 1, COL_REL_CN);

        if(hasRna)
        {
            addEntry(cells, widths, cellEntries, 1, COL_TPM);
            addEntry(cells, widths, cellEntries, 1, "Percentile");
            addEntry(cells, widths, cellEntries, 1, "Fold Change");
        }

        addEntry(cells, widths, cellEntries, 1, COL_DRIVER);

        Table table = Tables.createContent(width, intToFloatArray(widths), cellArray(cellEntries));

        for(PurpleGainDeletion gainDel : sort(gainsDels))
        {
            List<Cell> rowCells = Lists.newArrayList();

            rowCells.add(cells.createContent(gainDel.chromosome() + gainDel.chromosomeBand()));
            rowCells.add(cells.createContent(displayGene(gainDel)));
            rowCells.add(cells.createContent(gainDel.driver().type().toString()));
            rowCells.add(cells.createContent(gainDel.geneRange()));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(gainDel.minCopyNumber())));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(gainDel.maxCopyNumber())));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(gainDel.relativeCopyNumber())));

            if(hasRna)
            {
                rowCells.add(cells.createContent(formatTpmField(gainDel.tpm())));
                rowCells.add(cells.createContent(formatPercentileField(gainDel.tpmPercentile())));
                rowCells.add(cells.createContent(formatFoldChangeField(gainDel.tpmFoldChange())));
            }

            rowCells.add(cells.createContent(gainDel.driver().driverInterpretation().toString()));

            if(gainDel.driver().driverInterpretation() == DriverInterpretation.LOW)
            {
                reportResources.shadeCandidateCells(rowCells);
            }

            rowCells.forEach(x -> table.addCell(x));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    private static List<PurpleGainDeletion> sort(final List<PurpleGainDeletion> gainsAndDels)
    {
        return gainsAndDels.stream().sorted((gainDel1, gainDel2) ->
        {
            double likelihood1 = gainDel1.driver().driverLikelihood();
            double likelihood2 = gainDel2.driver().driverLikelihood();

            int likelihoodCompare = Double.compare(-likelihood1, -likelihood2);

            if(likelihoodCompare != 0)
                return likelihoodCompare;

            String location1 = Chromosomes.zeroPrefixed(gainDel1.chromosome() + gainDel1.chromosomeBand());
            String location2 = Chromosomes.zeroPrefixed(gainDel2.chromosome() + gainDel2.chromosomeBand());

            if(location1.equals(location2))
            {
                return gainDel1.gene().compareTo(gainDel2.gene());
            }
            else
            {
                return location1.compareTo(location2);
            }
        }).collect(Collectors.toList());
    }

    private static String displayGene(final PurpleGainDeletion gainDel)
    {
        String addon = Strings.EMPTY;
        if(!gainDel.isCanonical())
        {
            addon = " (alt)";
        }
        return gainDel.gene() + addon;
    }
}
