package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.algo.OrangeConstants.isCandidateLikelihood;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_CHR;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_CN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_REL_CN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_TYPE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.cellArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.intToFloatArray;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleChrArmCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;

public final class ChrArmCopyNumberTable
{
    public static Table build(
            final String title, float width, final List<PurpleChrArmCopyNumber> chrArmCopyNumbers, final ReportResources reportResources)
    {
        if(chrArmCopyNumbers.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<Cell> cellEntries = Lists.newArrayList();

        addEntry(cells, widths, cellEntries, 1, COL_CHR);
        addEntry(cells, widths, cellEntries, 1, "Arm");
        addEntry(cells, widths, cellEntries, 1, COL_TYPE);
        addEntry(cells, widths, cellEntries, 1, COL_CN);
        addEntry(cells, widths, cellEntries, 1, COL_REL_CN);
        addEntry(cells, widths, cellEntries, 1, COL_DRIVER);
        addEntry(cells, widths, cellEntries, 3, Strings.EMPTY); // to space things out

        Table table = Tables.createContent(width, intToFloatArray(widths), cellArray(cellEntries));

        for(PurpleChrArmCopyNumber chrArmCopyNumber : sort(chrArmCopyNumbers))
        {
            List<Cell> rowCells = Lists.newArrayList();

            rowCells.add(cells.createContent(chrArmCopyNumber.chromosome()));
            rowCells.add(cells.createContent(chrArmCopyNumber.arm()));
            rowCells.add(cells.createContent(chrArmCopyNumber.type()));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(chrArmCopyNumber.copyNumber())));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(chrArmCopyNumber.relativeCopyNumber())));
            rowCells.add(cells.createContent(chrArmCopyNumber.driverInterpretation().toString()));
            rowCells.add(cells.createContent(Strings.EMPTY));

            if(chrArmCopyNumber.driverInterpretation() == DriverInterpretation.LOW)
            {
                reportResources.shadeCandidateCells(rowCells);
            }

            rowCells.forEach(x -> table.addCell(x));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    private static List<PurpleChrArmCopyNumber> sort(final List<PurpleChrArmCopyNumber> arms)
    {
        return arms.stream().sorted((arm1, arm2) ->
        {
            DriverInterpretation driver1 = arm1.driverInterpretation();
            DriverInterpretation driver2 = arm2.driverInterpretation();

            if(driver1 != driver2)
            {
                return driver1 == DriverInterpretation.HIGH ? -1 : 1;
            }

            int chr1 = HumanChromosome.chromosomeRank(arm1.chromosome());
            int chr2 = HumanChromosome.chromosomeRank(arm2.chromosome());

            if(chr1 != chr2)
                return Integer.compare(chr1, chr2);

            return arm1.arm().compareTo(arm2.arm());

        }).collect(Collectors.toList());
    }

}
