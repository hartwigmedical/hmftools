package com.hartwig.hmftools.sage;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.sage.snpeff.SagePostProcess;
import com.hartwig.hmftools.sage.snpeff.SagePostProcessVCF;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class SagePostProcessApplication implements AutoCloseable {

    private static final Logger LOGGER = LogManager.getLogger(SagePostProcessApplication.class);

    private static final String IN_VCF = "in";
    private static final String OUT_VCF = "out";

    public static void main(String[] args) throws IOException, ParseException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String inputFilePath = cmd.getOptionValue(IN_VCF);
        final String outputFilePath = cmd.getOptionValue(OUT_VCF);

        if (outputFilePath == null || inputFilePath == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("SagePostProcessApplication", options);
            System.exit(1);
        }

        try (SagePostProcessApplication app = new SagePostProcessApplication(inputFilePath, outputFilePath)) {
            app.run();
        }
    }

    private final SagePostProcess postProcess;
    private final VCFFileReader fileReader;
    private final SagePostProcessVCF fileWriter;

    private SagePostProcessApplication(final String inputVCF, final String outputVCF) {
        LOGGER.info("Input: {}", inputVCF);
        LOGGER.info("Output: {}", outputVCF);

        this.fileReader = new VCFFileReader(new File(inputVCF), false);
        this.fileWriter = new SagePostProcessVCF(outputVCF);
        this.postProcess = new SagePostProcess(HmfGenePanelSupplier.allGenesMap37(), fileWriter::writeVariant);
    }

    private void run()  {
        fileWriter.writeHeader(fileReader.getFileHeader());

        for (VariantContext context : fileReader) {
            postProcess.accept(context);
        }
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(IN_VCF, true, "Input file.");
        options.addOption(OUT_VCF, true, "Output file.");
        return options;
    }

    @Override
    public void close() {
        fileReader.close();
        postProcess.close();
        fileWriter.close();

        LOGGER.info("Post process complete");
    }
}
