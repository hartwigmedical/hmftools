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

        String driverGenePanelFileRef = configBuilder.getValue(CFG_OLD_FILE_REF);
        String driverGenePanelFileNew = configBuilder.getValue(CFG_OLD_FILE_NEW);
        String outputDir = parseOutputDir(configBuilder);

        try
        {
            List<DriverGene> driverGenesRef = DriverGeneFile.read(driverGenePanelFileRef);
            List<DriverGene> driverGenesNew = DriverGeneFile.read(driverGenePanelFileNew);

            int refOnly = 0;
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

            for(DriverGene driverGeneRef : driverGenesRef)
            {
                DriverGene driverGeneNew = driverGenesNew.stream().filter(x -> x.gene().equals(driverGeneRef.gene())).findFirst().orElse(null);

                if(driverGeneNew == null)
                {
                    writeMissing(driverGeneRef, true, writer);
                    ++refOnly;
                }
                else
                {
                    matchedGenes.add(driverGeneRef.gene());

                    valueDiffs += checkDifferences(driverGeneRef, driverGeneNew, writer);
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

            CMP_LOGGER.info("additions: ref({}) new({}) valueDiffs({})", refOnly, newOnly, valueDiffs);
        }
        catch(Exception e)
        {
            CMP_LOGGER.error("failed to read driver gene panel file: {}", e.toString());
            System.exit(1);
        }
    }

    private static int checkDifferences(final DriverGene driverGeneRef, final DriverGene driverGeneNew, final BufferedWriter writer)
            throws IOException
    {
        String gene = driverGeneRef.gene();

        int diffs = 0;

        if(checkValueDifference(
                gene, "reportMissense", driverGeneRef.reportMissenseAndInframe(), driverGeneNew.reportMissenseAndInframe(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportNonsense", driverGeneRef.reportNonsenseAndFrameshift(), driverGeneNew.reportNonsenseAndFrameshift(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportSplice", driverGeneRef.reportSplice(), driverGeneNew.reportSplice(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportDeletion", driverGeneRef.reportDeletion(), driverGeneNew.reportDeletion(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportHetDeletion", driverGeneRef.reportHetDeletion(), driverGeneNew.reportHetDeletion(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportLoh", driverGeneRef.reportLoh(), driverGeneNew.reportLoh(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "hetDeletionThreshold", driverGeneRef.hetDeletionThreshold(), driverGeneNew.hetDeletionThreshold(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportDisruption", driverGeneRef.reportDisruption(), driverGeneNew.reportDisruption(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportAmplification", driverGeneRef.reportAmplification(), driverGeneNew.reportAmplification(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "amplificationRatio", driverGeneRef.amplificationRatio(), driverGeneNew.amplificationRatio(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportSomaticHotspot", driverGeneRef.reportSomaticHotspot(), driverGeneNew.reportSomaticHotspot(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "likelihoodType", driverGeneRef.likelihoodType().toString(), driverGeneNew.likelihoodType().toString(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportGermlineVariant", driverGeneRef.reportGermlineVariant().toString(),
                driverGeneNew.reportGermlineVariant().toString(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportGermlineHotspot", driverGeneRef.reportGermlineHotspot().toString(),
                driverGeneNew.reportGermlineHotspot().toString(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportGermlineDisruption", driverGeneRef.reportGermlineDisruption().toString(),
                driverGeneNew.reportGermlineDisruption().toString(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportGermlineAmplification", driverGeneRef.reportGermlineAmplification(), driverGeneNew.reportGermlineAmplification(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "additionalReportedTranscripts",
                driverGeneRef.additionalReportedTranscripts().stream().collect(Collectors.joining(ITEM_DELIM)),
                driverGeneNew.additionalReportedTranscripts().stream().collect(Collectors.joining(ITEM_DELIM)), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportPGX", driverGeneRef.reportPGX(), driverGeneNew.reportPGX(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportHighExpression", driverGeneRef.reportHighExpression(), driverGeneNew.reportHighExpression(), writer))
        {
            ++diffs;
        }

        if(checkValueDifference(
                gene, "reportLowExpression", driverGeneRef.reportLowExpression(), driverGeneNew.reportLowExpression(), writer))
        {
            ++diffs;
        }
        if(checkValueDifference(
                gene, "reportNovelSpliceJunction", driverGeneRef.reportNovelSpliceJunction(), driverGeneNew.reportNovelSpliceJunction(), writer))
        {
            ++diffs;
        }

        return diffs;
    }

    private static boolean checkValueDifference(
            final String gene, final String field, final boolean refValue, final boolean newValue, final BufferedWriter writer)
            throws IOException
    {
        if(refValue != newValue)
        {
            writeValueDifference(gene, field, String.valueOf(refValue), String.valueOf(newValue), writer);
            return true;
        }

        return false;
    }

    private static boolean checkValueDifference(
            final String gene, final String field, final double refValue, final double newValue, final BufferedWriter writer)
            throws IOException
    {
        if(refValue != newValue)
        {
            writeValueDifference(gene, field, String.valueOf(refValue), String.valueOf(newValue), writer);
            return true;
        }

        return false;
    }

    private static boolean checkValueDifference(
            final String gene, final String field, final String refValue, final String newValue, final BufferedWriter writer)
            throws IOException
    {
        if(!refValue.equals(newValue))
        {
            writeValueDifference(gene, field, refValue, newValue, writer);
            return true;
        }

        return false;
    }

    private static void writeValueDifference(final String gene, final String field, final String refValue, final String newValue, final BufferedWriter writer)
            throws IOException
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(gene);
        sj.add("VALUE");
        sj.add(format("%s(%s/%s)", field, refValue, newValue));
        writer.write(sj.toString());
        writer.newLine();
    }

    private static void writeMissing(final DriverGene driverGene, boolean isRef, final BufferedWriter writer)
            throws IOException
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(driverGene.gene());
        sj.add(isRef ? "REF_ONLY" : "NEW_ONLY");
        sj.add("");
        writer.write(sj.toString());
        writer.newLine();
    }

    private static final String CFG_OLD_FILE_REF = "driver_gene_panel_ref";
    private static final String CFG_OLD_FILE_NEW = "driver_gene_panel_new";

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder("Compar");

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
        configBuilder.addPath(CFG_OLD_FILE_REF, true, "Path to reference driver gene panel");
        configBuilder.addPath(CFG_OLD_FILE_NEW, true, "Path to new driver gene panel");

        configBuilder.checkAndParseCommandLine(args);

        new ComparDriverGeneFiles(configBuilder);
    }
}
