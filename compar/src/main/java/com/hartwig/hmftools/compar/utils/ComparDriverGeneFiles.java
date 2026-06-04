package com.hartwig.hmftools.compar.utils;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class ComparDriverGeneFiles
{
    public ComparDriverGeneFiles(final ConfigBuilder configBuilder)
    {
        CMP_LOGGER.info("comparing driver gene panel files");

        String driverGenePanelFileOld = configBuilder.getValue(CFG_PANEL_FILE_OLD);
        String driverGenePanelFileNew = configBuilder.getValue(CFG_PANEL_FILE_NEW);
        String outputDir = parseOutputDir(configBuilder);

        try
        {
            List<DriverGene> driverGenesOld = DriverGeneFile.read(driverGenePanelFileOld);
            List<DriverGene> driverGenesNew = DriverGeneFile.read(driverGenePanelFileNew);

            int oldOnly = 0;
            int newOnly = 0;
            int valueDiffs = 0;

            String outputFile = outputDir + "compar_driver_gene_panels.tsv";
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(FLD_GENE_NAME);
            sj.add("MismatchType");
            sj.add("Difference");

            writer.write(sj.toString());
            writer.newLine();

            Set<String> matchedGenes = Sets.newHashSet();

            for(DriverGene driverGeneOld : driverGenesOld)
            {
                DriverGene driverGeneNew = driverGenesNew.stream().filter(x -> x.gene().equals(driverGeneOld.gene())).findFirst().orElse(null);

                if(driverGeneNew == null)
                {
                    writeMissing(driverGeneOld, true, writer);
                    ++oldOnly;
                }
                else
                {
                    matchedGenes.add(driverGeneOld.gene());

                    valueDiffs += checkDifferences(driverGeneOld, driverGeneNew, writer);
                }
            }

            for(DriverGene driverGeneNew : driverGenesNew)
            {
                if(!matchedGenes.contains(driverGeneNew.gene()))
                {
                    writeMissing(driverGeneNew, false, writer);
                    ++newOnly;
                }
            }

            writer.close();

            CMP_LOGGER.info("additions: old({}) new({}) valueDiffs({})", oldOnly, newOnly, valueDiffs);
        }
        catch(Exception e)
        {
            CMP_LOGGER.error("failed to read driver gene panel file: {}", e.toString());
            System.exit(1);
        }
    }

    private static int checkDifferences(final DriverGene driverGeneOld, final DriverGene driverGeneNew, final BufferedWriter writer)
            throws IOException
    {
        String gene = driverGeneOld.gene();

        int diffs = 0;

        if(checkValueDifference(
                gene, "reportMissense", driverGeneOld.reportMissenseAndInframe(), driverGeneNew.reportMissenseAndInframe(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportNonsense", driverGeneOld.reportNonsenseAndFrameshift(), driverGeneNew.reportNonsenseAndFrameshift(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportSplice", driverGeneOld.reportSplice(), driverGeneNew.reportSplice(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportDeletion", driverGeneOld.reportDeletion(), driverGeneNew.reportDeletion(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportHetDeletion", driverGeneOld.reportHetDeletion(), driverGeneNew.reportHetDeletion(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportLoh", driverGeneOld.reportLoh(), driverGeneNew.reportLoh(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "hetDeletionThreshold", driverGeneOld.hetDeletionThreshold(), driverGeneNew.hetDeletionThreshold(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportDisruption", driverGeneOld.reportDisruption(), driverGeneNew.reportDisruption(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportAmplification", driverGeneOld.reportAmplification(), driverGeneNew.reportAmplification(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "amplificationRatio", driverGeneOld.amplificationRatio(), driverGeneNew.amplificationRatio(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportSomaticHotspot", driverGeneOld.reportSomaticHotspot(), driverGeneNew.reportSomaticHotspot(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "likelihoodType", driverGeneOld.likelihoodType().toString(), driverGeneNew.likelihoodType().toString(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportGermlineVariant", driverGeneOld.reportGermlineVariant().toString(),
                driverGeneNew.reportGermlineVariant().toString(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportGermlineHotspot", driverGeneOld.reportGermlineHotspot().toString(),
                driverGeneNew.reportGermlineHotspot().toString(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportGermlineDisruption", driverGeneOld.reportGermlineDisruption().toString(),
                driverGeneNew.reportGermlineDisruption().toString(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportGermlineAmplification", driverGeneOld.reportGermlineAmplification(), driverGeneNew.reportGermlineAmplification(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "additionalReportedTranscripts",
                driverGeneOld.additionalReportedTranscripts().stream().collect(Collectors.joining(ITEM_DELIM)),
                driverGeneNew.additionalReportedTranscripts().stream().collect(Collectors.joining(ITEM_DELIM)), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportPGX", driverGeneOld.reportPGX(), driverGeneNew.reportPGX(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportHighExpression", driverGeneOld.reportHighExpression(), driverGeneNew.reportHighExpression(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportLowExpression", driverGeneOld.reportLowExpression(), driverGeneNew.reportLowExpression(), writer))
        {
            ++diffs;
        }
        if(checkValueDifference(
                gene, "reportNovelSpliceJunction", driverGeneOld.reportNovelSpliceJunction(), driverGeneNew.reportNovelSpliceJunction(), writer))
        {
            ++diffs;
        }

        return diffs;
    }

    private static boolean checkValueDifference(
            final String gene, final String field, final boolean oldValue, final boolean newValue, final BufferedWriter writer)
            throws IOException
    {
        if(oldValue != newValue)
        {
            writeValueDifference(gene, field, String.valueOf(oldValue), String.valueOf(newValue), writer);
            return true;
        }

        return false;
    }

    private static boolean checkValueDifference(
            final String gene, final String field, final double oldValue, final double newValue, final BufferedWriter writer)
            throws IOException
    {
        if(oldValue != newValue)
        {
            writeValueDifference(gene, field, String.valueOf(oldValue), String.valueOf(newValue), writer);
            return true;
        }

        return false;
    }

    private static boolean checkValueDifference(
            final String gene, final String field, final String oldValue, final String newValue, final BufferedWriter writer)
            throws IOException
    {
        if(!oldValue.equals(newValue))
        {
            writeValueDifference(gene, field, oldValue, newValue, writer);
            return true;
        }

        return false;
    }

    private static void writeValueDifference(final String gene, final String field, final String oldValue, final String newValue, final BufferedWriter writer)
            throws IOException
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(gene);
        sj.add("VALUE");
        sj.add(format("%s(%s/%s)", field, oldValue, newValue));
        writer.write(sj.toString());
        writer.newLine();
    }

    private static void writeMissing(final DriverGene driverGene, boolean isOld, final BufferedWriter writer)
            throws IOException
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(driverGene.gene());
        sj.add(isOld ? "OLD_ONLY" : "NEW_ONLY");
        sj.add("");
        writer.write(sj.toString());
        writer.newLine();
    }

    private static final String CFG_PANEL_FILE_OLD = "driver_gene_panel_old";
    private static final String CFG_PANEL_FILE_NEW = "driver_gene_panel_new";

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder("Compar");

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
        configBuilder.addPath(CFG_PANEL_FILE_OLD, true, "Path to old driver gene panel");
        configBuilder.addPath(CFG_PANEL_FILE_NEW, true, "Path to new driver gene panel");

        configBuilder.checkAndParseCommandLine(args);

        new ComparDriverGeneFiles(configBuilder);
    }
}
