package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.formatPercentileField;
import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.ReportResources.formatTpmField;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.datamodel.isofox.IsofoxRecord;
import com.hartwig.hmftools.datamodel.isofox.NovelSpliceJunction;
import com.hartwig.hmftools.datamodel.isofox.RnaFusion;
import com.hartwig.hmftools.datamodel.linx.FusionLikelihoodType;
import com.hartwig.hmftools.datamodel.linx.FusionPhasedType;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.linx.LinxFusionType;
import com.hartwig.hmftools.datamodel.linx.LinxUnreportableReason;
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
        Table table = Tables.createContent(width,
                new float[] { 1, 5 },
                new Cell[] { cells.createHeader("Fusion"), cells.createHeader("Details") });

        for(LinxFusion fusion : sortLinxFusions(fusions))
        {
            table.addCell(cells.createContent(fusion.display()));

            Table details = new Table(UnitValue.createPercentArray(new float[] { 1, 3 }));
            Stream.of(Maps.immutableEntry("5' End", cells.createValue(fiveEndString(fusion))),
                            Maps.immutableEntry("3' Start", cells.createValue(threeStartString(fusion))),
                            Maps.immutableEntry("Junction CN", cells.createValue(formatSingleDigitDecimal(fusion.junctionCopyNumber()))),
                            Maps.immutableEntry("RNA support",
                                    cells.createValue(rnaFragmentSupportTable(isofox, fusion, cells)).setKeepTogether(true)),
                            Maps.immutableEntry("Phasing", cells.createValue(display(fusion.phased()))),
                            Maps.immutableEntry("Reported type (DL)",
                                    cells.createValue(fusion.reportedType() + " (" + display(fusion.driverLikelihood()) + ")")),
                            Maps.immutableEntry("Unreported reason(s)", cells.createValue(display(fusion.unreportedReasons()))),
                            Maps.immutableEntry("Chain links (terminated?)",
                                    cells.createValue(fusion.chainLinks() + (fusion.chainTerminated() ? " (Yes)" : " (No)"))),
                            Maps.immutableEntry("Domains kept", cells.createValue(!fusion.domainsKept().isEmpty() ? fusion.domainsKept() : "-")),
                            Maps.immutableEntry("Domains lost", cells.createValue(!fusion.domainsLost().isEmpty() ? fusion.domainsLost() : "-")))
                    .forEach(entry ->
                    {
                        details.addCell(cells.createKey(entry.getKey()));
                        details.addCell(entry.getValue());
                    });
            // Need to keep this details table to avoid page-wrapping that cuts through the middle of a single fusion
            table.addCell(cells.createContent(details).setKeepTogether(true));
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

    private static String display(FusionLikelihoodType fusionLikelihoodType)
    {
        switch(fusionLikelihoodType)
        {
            case HIGH:
                return "High";
            case LOW:
                return "Low";
            case NA:
                return "NA";
        }
        throw new IllegalStateException();
    }

    private static String display(final List<LinxUnreportableReason> unreportedReasons)
    {
        return unreportedReasons.stream().map(item ->
        {
            switch(item)
            {
                case NONE:
                    return "-";
                case NOT_KNOWN:
                    return "Not a known fusion pair";
                case UNPHASED_NOT_KNOWN:
                    return "Unphased, not a known fusion pair";
                case UNPHASED_5P_UTR:
                    return "Unphased, 5' UTR";
                case UNPHASED_SHORT:
                    return "Unphased, short unphased distance";
                case SGL_NOT_KNOWN:
                    return "SGL, no known fusion pair";
                case PRE_GENE_DISTANCE:
                    return "Max upstream distance exceeded";
                case NONSENSE_MEDIATED_DECAY:
                    return "Nonsense mediated decay";
                case NEG_SPLICE_ACC_DISTANCE:
                    return "Negative previous splice acceptor distance";
                case EXON_SKIPPING:
                    return "Exon skipping";
                case CHAIN_TERMINATED:
                    return "Chain terminated";
                case NON_DISRUPTIVE_CHAIN:
                    return "Non-disruptive chain";
                case INVALID_TRAVERSAL:
                    return "Invalid chain traversal";
                case CHAIN_LINKS:
                    return "Maximum chain links exceeded";
                case DISRUPTED_PROTEIN_DOMAINS:
                    return "Disrupted protein domains";
                default:
                    throw new IllegalArgumentException("Unknown unreportable reason: " + item);
            }
        }).collect(Collectors.joining(", "));
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
        // String fcDbString = "FC " + Expressions.foldChangeDatabase(geneEndExpression);
        // String dbString = " DB percentile " + Expressions.percentileDatabase(geneEndExpression) + " (" + fcDbString + ")";

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
            if(fusion1.driverLikelihood() == fusion2.driverLikelihood())
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
                return fusion1.driverLikelihood() == FusionLikelihoodType.HIGH ? -1 : 1;
            }
        }).collect(Collectors.toList());
    }

    private static <T> List<T> max5(final List<T> elements)
    {
        return elements.subList(0, Math.min(5, elements.size()));
    }
}
