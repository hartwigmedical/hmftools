package com.hartwig.hmftools.sage;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.CanonicalTranscript;
import com.hartwig.hmftools.common.genome.region.CanonicalTranscriptFactory;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.sage.SagePostProcessVCF;
import com.hartwig.hmftools.sage.snpeff.SagePostProcess;

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

    private static final String HG19 = "hg19";
    private static final String HG38 = "hg38";

    private static final Logger LOGGER = LogManager.getLogger(SagePostProcessApplication.class);

    private static final String IN_VCF = "in";
    private static final String OUT_VCF = "out";
    private static final String ASSEMBLY = "assembly";

    public static void main(String[] args) throws IOException, ParseException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String inputFilePath = cmd.getOptionValue(IN_VCF);
        final String outputFilePath = cmd.getOptionValue(OUT_VCF);
        final String assembly = cmd.getOptionValue(ASSEMBLY);

        if (assembly == null || !assembly.equals(HG19) && !assembly.equals(HG38)) {
            LOGGER.error("Parameter assembly must be one of [hg19, hg38]");
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("SagePostProcessApplication", options);
            System.exit(1);
        }

        if (outputFilePath == null || inputFilePath == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("SagePostProcessApplication", options);
            System.exit(1);
        }

        try (SagePostProcessApplication app = new SagePostProcessApplication(inputFilePath, outputFilePath, assembly)) {
            app.run();
        }
    }

    private final SagePostProcess postProcess;
    private final VCFFileReader fileReader;
    private final SagePostProcessVCF fileWriter;

    private SagePostProcessApplication(final String inputVCF, final String outputVCF, final String assembly) {
        LOGGER.info("Input: {}", inputVCF);
        LOGGER.info("Output: {}", outputVCF);

        this.fileReader = new VCFFileReader(new File(inputVCF), false);
        this.fileWriter = new SagePostProcessVCF(outputVCF);
        Map<String, HmfTranscriptRegion> geneMap =
                assembly.equals(HG19) ? HmfGenePanelSupplier.allGenesMap37() : HmfGenePanelSupplier.allGenesMap38();
        List<CanonicalTranscript> canonicalTranscripts =
                assembly.equals(HG19) ? CanonicalTranscriptFactory.create37() : CanonicalTranscriptFactory.create38();
        this.postProcess = new SagePostProcess(geneMap, canonicalTranscripts, fileWriter::writeVariant);
    }

    private void run() {
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
        options.addOption(ASSEMBLY, true, "Assembly. Must be one of [hg19, hg38]");
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
