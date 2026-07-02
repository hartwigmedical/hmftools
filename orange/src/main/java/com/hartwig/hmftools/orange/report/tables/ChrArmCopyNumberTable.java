package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_CHR;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_CN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_REL_CN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_TYPE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.stringArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createStandardTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createEmptyTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.toPercentages;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.intToFloatArray;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleChrArmCopyNumber;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;

import be.quodlibet.boxable.Cell;
import be.quodlibet.boxable.BaseTable;

import org.apache.pdfbox.pdmodel.PDPage;

import com.hartwig.hmftools.orange.report.DocumentContext;

import java.io.IOException;

import org.apache.logging.log4j.util.Strings;

public final class ChrArmCopyNumberTable
{
    public static BaseTable build(final DocumentContext docCtx,
            final String title, float width, final List<PurpleChrArmCopyNumber> chrArmCopyNumbers, final ReportResources reportResources)
            throws IOException
    {
        if(chrArmCopyNumbers.isEmpty())
        {
            return createEmptyTable(docCtx, title, width, reportResources);
        }

        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<String> headers = Lists.newArrayList();

        addEntry(widths, headers, 1, COL_CHR);
        addEntry(widths, headers, 1, "Arm");
        addEntry(widths, headers, 1, COL_TYPE);
        addEntry(widths, headers, 1, COL_CN);
        addEntry(widths, headers, 1, COL_REL_CN);
        addEntry(widths, headers, 1, COL_DRIVER);
        addEntry(widths, headers, 3, Strings.EMPTY); // to space things out

        BaseTable table = createStandardTable(docCtx, title, width, intToFloatArray(widths), stringArray(headers), reportResources);
        float[] pcts = toPercentages(intToFloatArray(widths));

        for(PurpleChrArmCopyNumber chrArmCopyNumber : sort(chrArmCopyNumbers))
        {
            List<String> rowCells = Lists.newArrayList();

            rowCells.add(cells.createContent(chrArmCopyNumber.chromosome()));
            rowCells.add(cells.createContent(chrArmCopyNumber.arm()));
            rowCells.add(cells.createContent(chrArmCopyNumber.type()));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(chrArmCopyNumber.copyNumber())));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(chrArmCopyNumber.relativeCopyNumber())));
            rowCells.add(cells.createContent(chrArmCopyNumber.driverInterpretation().toString()));
            rowCells.add(cells.createContent(Strings.EMPTY));

            List<Cell<PDPage>> createdCells = cells.addRow(table, pcts, rowCells);

            if(chrArmCopyNumber.driverInterpretation() == DriverInterpretation.LOW)
            {
                reportResources.shadeCandidateCells(createdCells);
            }
        }

        return table;
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
            {
                return Integer.compare(chr1, chr2);
            }

            return arm1.arm().compareTo(arm2.arm());

        }).collect(Collectors.toList());
    }
}
