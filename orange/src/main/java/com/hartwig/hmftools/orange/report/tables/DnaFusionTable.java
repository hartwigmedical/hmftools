package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.formatPercentileField;
import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.ReportResources.formatTpmField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_JCN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_RNA;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.cellArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.floatArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.intToFloatArray;

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
        addEntry(cells, widths, cellEntries, 2, "5' End");
        addEntry(cells, widths, cellEntries, 2, "3' End");
        addEntry(cells, widths, cellEntries, 0.8, COL_JCN);
        addEntry(cells, widths, cellEntries, 1.2, "Phasing");
        addEntry(cells, widths, cellEntries, 2, "Type");
        addEntry(cells, widths, cellEntries, 2, "Domains kept");
        addEntry(cells, widths, cellEntries, 2, "Domains lost");

        if(hasRna)
            addEntry(cells, widths, cellEntries, 1, COL_RNA);

        addEntry(cells, widths, cellEntries, 1, COL_DRIVER);

        Table table = Tables.createContent(width, floatArray(widths), cellArray(cellEntries));

        for(LinxFusion fusion : sortLinxFusions(fusions))
        {
            List<Cell> rowCells = Lists.newArrayList();

            rowCells.add(cells.createContent(fusion.display()));
            rowCells.add(cells.createValue(fiveEndString(fusion)));
            rowCells.add(cells.createValue(threeStartString(fusion)));
            rowCells.add(cells.createValue(formatSingleDigitDecimal(fusion.junctionCopyNumber())));
            rowCells.add(cells.createValue(display(fusion.phased())));
            rowCells.add(cells.createValue(fusion.reportedType().toString()));
            rowCells.add(cells.createValue(!fusion.domainsKept().isEmpty() ? fusion.domainsKept() : "-"));
            rowCells.add(cells.createValue(!fusion.domainsLost().isEmpty() ? fusion.domainsLost() : "-"));

            if(hasRna)
                rowCells.add(cells.createValue(rnaFragmentSupportTable(isofox, fusion, cells)));

            rowCells.add(cells.createValue(fusion.driverInterpretation().toString()));

            if(fusion.driverInterpretation() == DriverInterpretation.LOW)
            {
                reportResources.shadeCandidateCells(rowCells);
            }

            rowCells.forEach(x -> table.addCell(x));
        }

        return new Tables(reportResources).createWrapping(table, title);
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

    private static String fiveEndString(final LinxFusion fusion)
    {
        return fusion.geneStart() + " " + fusion.geneContextStart() + " (" + fusion.geneTranscriptStart() + ")";
    }

    private static String threeStartString(final LinxFusion fusion)
    {
        return fusion.geneEnd() + " " + fusion.geneContextEnd() + " (" + fusion.geneTranscriptEnd() + ")";
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
        else if(fusion.geneStart().equals(fusion.geneEnd()))
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
        GeneExpression geneEndExpression = Expressions.findByGene(isofox.allGeneExpressions(), fusion.geneEnd());

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
            if(junction.gene().equals(fusion.geneStart()) && junction.gene().equals(fusion.geneEnd()))
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
                if(fusion1.geneStart().equals(fusion2.geneStart()))
                {
                    return fusion1.geneEnd().compareTo(fusion2.geneEnd());
                }
                else
                {
                    return fusion1.geneStart().compareTo(fusion2.geneStart());
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
