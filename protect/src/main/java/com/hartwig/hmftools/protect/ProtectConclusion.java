package com.hartwig.hmftools.protect;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.Variant;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;
import com.hartwig.hmftools.protect.common.GenomicData;
import com.hartwig.hmftools.protect.report.chord.ChordAnalysis;
import com.hartwig.hmftools.protect.report.chord.ChordFileReader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ProtectConclusion {
    private static final Logger LOGGER = LogManager.getLogger(ProtectActionability.class);

    private static final String TUMOR_SAMPLE_ID = "tumor_sample_id";
    private static final String SOMATIC_VARIANT_VCF = "somatic_variant_vcf";
    private static final String PURPLE_PURITY_TSV = "purple_purity_tsv";
    private static final String PURPLE_GENE_CNV_TSV = "purple_gene_cnv_tsv";
    private static final String LINX_FUSION_TSV = "linx_fusion_tsv";
    private static final String CHORD_TXT = "chord_txt";


    private static final String CONCLUSION_TSV = "conclusion_tsv";

    public static void main(@NotNull final String[] args) throws ParseException, IOException {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        String tumorSampleId = cmd.getOptionValue(TUMOR_SAMPLE_ID);

        // Params specific for specific sample
        final String somaticVariantVcf = cmd.getOptionValue(SOMATIC_VARIANT_VCF);
        final String purplePurityTsv = cmd.getOptionValue(PURPLE_PURITY_TSV);
        final String purpleGeneCnvTsv = cmd.getOptionValue(PURPLE_GENE_CNV_TSV);
        final String linxFusionTsv = cmd.getOptionValue(LINX_FUSION_TSV);
        final String chordTxt = cmd.getOptionValue(CHORD_TXT);


        final String OutputConclusionTsv = cmd.getOptionValue(CONCLUSION_TSV);

        if (!validInputForBaseReport(cmd)) {
            printUsageAndExit(options);
        }

        List<? extends Variant> passSomaticVariants = GenomicData.readPassSomaticVariants(tumorSampleId, somaticVariantVcf);
        double ploidy = GenomicData.extractPloidy(purplePurityTsv);
        List<GeneCopyNumber> geneCopyNumbers = GenomicData.readGeneCopyNumbers(purpleGeneCnvTsv);
        List<ReportableGeneFusion> geneFusions = GenomicData.readGeneFusions(linxFusionTsv);

        int tumorMTL = 0;
        int tumorMSI = 0;
        double chordScore = ChordFileReader.read(chordTxt).hrdValue();
        int tumorMTB = 0;

        LOGGER.info("Create conclusion for sample");
        writeConclusionOfSample(OutputConclusionTsv, "");

    }

    private static void writeConclusionOfSample(@NotNull String OutputConclusionTsv, @NotNull String conclusion) throws IOException {
        //TODO create conclusion
        BufferedWriter writerReport = new BufferedWriter(new FileWriter(OutputConclusionTsv, false));

        writerReport.write(conclusion);

        writerReport.close();
    }

    private static boolean validInputForBaseReport(@NotNull CommandLine cmd) {
        return valueExists(cmd, TUMOR_SAMPLE_ID) && fileExists(cmd, SOMATIC_VARIANT_VCF) && fileExists(cmd, PURPLE_PURITY_TSV)
                && fileExists(cmd, PURPLE_GENE_CNV_TSV) && fileExists(cmd, LINX_FUSION_TSV) && fileExists(cmd, CHORD_TXT);
    }

    private static boolean valueExists(@NotNull CommandLine cmd, @NotNull String param) {
        String value = cmd.getOptionValue(param);
        if (value == null) {
            LOGGER.warn(param + " has to be provided");
            return false;
        }
        return true;
    }

    private static boolean fileExists(@NotNull CommandLine cmd, @NotNull String param) {
        String value = cmd.getOptionValue(param);

        if (value == null || !pathExists(value)) {
            LOGGER.warn(param + " has to be an existing file: " + value);
            return false;
        }

        return true;
    }

    private static boolean pathExists(@NotNull String path) {
        return Files.exists(new File(path).toPath());
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();

        options.addOption(TUMOR_SAMPLE_ID, true, "The sample ID for which a patient report will be generated.");

        options.addOption(SOMATIC_VARIANT_VCF, true, "Path towards the somatic variant VCF.");
        options.addOption(PURPLE_PURITY_TSV, true, "Path towards the purple purity TSV.");
        options.addOption(PURPLE_GENE_CNV_TSV, true, "Path towards the purple gene copy number TSV.");
        options.addOption(LINX_FUSION_TSV, true, "Path towards the linx fusion TSV.");
        options.addOption(CHORD_TXT, true, "Path towards the chord txt file.");

        options.addOption(CONCLUSION_TSV, true, "Path towards the conclusion TSV.");

        return parser.parse(options, args);
    }

    @NotNull
    private static Options createBasicOptions() {
        final Options options = new Options();

        return options;
    }

    private static void printUsageAndExit(@NotNull Options options) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Protect Conclusion", options);
        System.exit(1);
    }

}
