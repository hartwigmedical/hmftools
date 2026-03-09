package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import static com.hartwig.hmftools.orange.report.ReportResources.formatPercentileField;
import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.ReportResources.formatTpmField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_JCN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_RNA;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.cellArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.floatArray;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.datamodel.isofox.IsofoxRecord;
import com.hartwig.hmftools.datamodel.isofox.NovelSpliceJunction;
import com.hartwig.hmftools.datamodel.isofox.RnaFusion;
import com.hartwig.hmftools.datamodel.linx.FusionPhasedType;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.linx.LinxFusionType;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Expressions;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.IBlockElement;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.Nullable;

public final class DnaFusionTable
{
    public static Table build(
            final String title, float width, final List<LinxFusion> fusions, @Nullable final IsofoxRecord isofox,
            final ReportResources reportResources)
    {
        if(fusions.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);

        List<Float> widths = Lists.newArrayList();
        List<Cell> cellEntries = Lists.newArrayList();

        boolean hasRna = isofox != null;

        addEntry(cells, widths, cellEntries, 2, "Fusion");
        addEntry(cells, widths, cellEntries, 5, "Junctions");
        addEntry(cells, widths, cellEntries, 1, COL_JCN);
        addEntry(cells, widths, cellEntries, 1, "Phasing");
        addEntry(cells, widths, cellEntries, 2, "Type");

        if(hasRna)
            addEntry(cells, widths, cellEntries, 2, COL_RNA);

        addEntry(cells, widths, cellEntries, 1, COL_DRIVER);

        Table table = Tables.createContent(width, floatArray(widths), cellArray(cellEntries));

        for(LinxFusion fusion : sortLinxFusions(fusions))
        {
            List<Cell> rowCells = Lists.newArrayList();

            rowCells.add(cells.createContent(fusionDisplay(fusion)));
            rowCells.add(cells.createContent(transcriptJunctions(fusion)));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(fusion.junctionCopyNumber())));
            rowCells.add(cells.createContent(display(fusion.phased())));
            rowCells.add(cells.createContent(fusion.reportedType().toString()));

            /*
            if(hasRna)
                rowCells.add(cells.createContent(rnaFragmentSupportTable(isofox, fusion, cells)));
            */

            rowCells.add(cells.createContent(format("SR=%d DF=%d", 10, 30)));

            rowCells.add(cells.createContent(fusion.driverInterpretation().toString()));

            if(fusion.driverInterpretation() == DriverInterpretation.LOW)
            {
                reportResources.shadeCandidateCells(rowCells);
            }

            rowCells.forEach(x -> table.addCell(x));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    private static String fusionDisplay(final LinxFusion fusion)
    {
        return format("%s::%s", fusion.geneUp(), fusion.geneDown());
    }

    private static String transcriptJunctions(final LinxFusion fusion)
    {
        return format("%s (%s) - %s (%s)",
                fusion.contextUp(), fusion.transcriptUp(), fusion.contextDown(), fusion.transcriptDown());

    }

    private static String display(FusionPhasedType fusionPhasedType)
    {
        switch(fusionPhasedType)
        {
            case INFRAME:
                return "Inframe";
            case SKIPPED_EXONS:
                return "Skipped exons";
            case OUT_OF_FRAME:
                return "Out of frame";
        }
        throw new IllegalStateException();
    }

    private static IBlockElement rnaFragmentSupportTable(@Nullable final IsofoxRecord isofox, final LinxFusion fusion, final Cells cells)
    {
        if(isofox == null)
        {
            return new Paragraph(ReportResources.NOT_AVAILABLE);
        }

        if(fusion.reportedType().equals(LinxFusionType.IG_KNOWN_PAIR) || fusion.reportedType().equals(LinxFusionType.IG_PROMISCUOUS))
        {
            return supportFromExpressionOfGeneEnd(isofox, fusion);
        }
        else if(fusion.geneUp().equals(fusion.geneDown()))
        {
            return supportFromSpliceJunctions(isofox, fusion, cells);
        }
        else
        {
            return supportFromRnaFusions(isofox, fusion, cells);
        }
    }

    private static IBlockElement supportFromExpressionOfGeneEnd(final IsofoxRecord isofox, final LinxFusion fusion)
    {
        GeneExpression geneEndExpression = Expressions.findByGene(isofox.allGeneExpressions(), fusion.geneDown());

        if(geneEndExpression == null)
        {
            return new Paragraph("None");
        }

        String tpmString = "TPM " + formatTpmField(geneEndExpression.tpm());
        String fcTypeString = "FC " + Expressions.formatFoldChange(geneEndExpression);
        String typeString = " Type percentile " + formatPercentileField(geneEndExpression.percentileCohort()) + " (" + fcTypeString + ")";

        return new Paragraph(geneEndExpression.gene() + " " + tpmString + ", " + typeString);
    }

    private static IBlockElement supportFromSpliceJunctions(final IsofoxRecord isofox, final LinxFusion fusion, final Cells cells)
    {
        List<NovelSpliceJunction> matches = Lists.newArrayList();
        for(NovelSpliceJunction junction : isofox.allNovelSpliceJunctions())
        {
            if(junction.gene().equals(fusion.geneUp()) && junction.gene().equals(fusion.geneDown()))
            {
                matches.add(junction);
            }
        }

        if(matches.isEmpty())
        {
            return new Paragraph("None");
        }

        Table fragmentSupportTable = new Table(UnitValue.createPercentArray(new float[] { 1 }));
        for(NovelSpliceJunction junction : max5(sortNovelSpliceJunctions(matches)))
        {
            String position = junction.chromosome() + ":" + junction.junctionStart() + "-" + junction.junctionEnd();
            String fragments = junction.fragmentCount() + " fragments";
            String depth = junction.depthStart() + " / " + junction.depthEnd() + " depth";
            fragmentSupportTable.addCell(cells.createValue(position + ", " + junction.type() + " (" + fragments + ", " + depth + ")"));
        }

        return fragmentSupportTable;
    }

    private static List<NovelSpliceJunction> sortNovelSpliceJunctions(final List<NovelSpliceJunction> novelSpliceJunctions)
    {
        return novelSpliceJunctions.stream()
                .sorted((junction1, junction2) -> Integer.compare(junction2.fragmentCount(), junction1.fragmentCount()))
                .collect(Collectors.toList());
    }

    private static IBlockElement supportFromRnaFusions(final IsofoxRecord isofox, final LinxFusion fusion, final Cells cells)
    {
        List<RnaFusion> matches = Lists.newArrayList();
        for(RnaFusion rnaFusion : isofox.allFusions())
        {
            if(rnaFusion.display().equals(fusion.display()))
            {
                matches.add(rnaFusion);
            }
        }

        if(matches.isEmpty())
        {
            return new Paragraph("None");
        }

        Table fragmentSupportTable = new Table(UnitValue.createPercentArray(new float[] { 1 }));
        for(RnaFusion rnaFusion : max5(sortRnaFusions(matches)))
        {
            String up = rnaFusion.chromosomeStart() + ":" + rnaFusion.positionStart();
            String down = rnaFusion.chromosomeEnd() + ":" + rnaFusion.positionEnd();
            String position = up + "-" + down;

            String split = rnaFusion.splitFragments() + " split";
            String realigned = rnaFusion.realignedFrags() + " realig.";
            String discord = rnaFusion.discordantFrags() + " discord.";
            String fragments = split + " / " + realigned + " / " + discord + " fragments";

            String depth = rnaFusion.depthStart() + " / " + rnaFusion.depthEnd() + " depth";
            fragmentSupportTable.addCell(cells.createValue(position + ", " + fragments + ", " + depth));
        }

        return fragmentSupportTable;
    }

    private static List<RnaFusion> sortRnaFusions(final List<RnaFusion> rnaFusions)
    {
        return rnaFusions.stream().sorted((fusion1, fusion2) ->
        {
            int sumFragments1 = fusion1.splitFragments() + fusion1.realignedFrags() + fusion1.discordantFrags();
            int sumFragments2 = fusion2.splitFragments() + fusion2.realignedFrags() + fusion2.discordantFrags();
            return Integer.compare(sumFragments2, sumFragments1);
        }).collect(Collectors.toList());
    }

    private static List<LinxFusion> sortLinxFusions(final List<LinxFusion> fusions)
    {
        return fusions.stream().sorted((fusion1, fusion2) ->
        {
            if(fusion1.driverInterpretation() == fusion2.driverInterpretation())
            {
                if(fusion1.geneUp().equals(fusion2.geneUp()))
                {
                    return fusion1.geneDown().compareTo(fusion2.geneDown());
                }
                else
                {
                    return fusion1.geneUp().compareTo(fusion2.geneUp());
                }
            }
            else
            {
                return fusion1.driverInterpretation() == DriverInterpretation.HIGH ? -1 : 1;
            }
        }).collect(Collectors.toList());
    }

    private static <T> List<T> max5(final List<T> elements)
    {
        return elements.subList(0, Math.min(5, elements.size()));
    }
}
