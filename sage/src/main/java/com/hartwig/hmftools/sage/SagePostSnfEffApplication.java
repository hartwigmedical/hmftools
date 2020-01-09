package com.hartwig.hmftools.sage;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.sage.snpeff.Reannotate;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class SagePostSnfEffApplication implements AutoCloseable {

    private static final String IN_VCF = "in";
    private static final String OUT_VCF = "out";

    public static void main(String[] args) throws IOException, ParseException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String inputFilePath = cmd.getOptionValue(IN_VCF);
        final String outputFilePath = cmd.getOptionValue(OUT_VCF);

        if (outputFilePath == null || inputFilePath == null ) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("SageHotspotAnnotation", options);
            System.exit(1);
        }

        try (SagePostSnfEffApplication app = new SagePostSnfEffApplication(inputFilePath)) {
            app.run();
        }
    }

    private final String inputVcf;
    private final Reannotate reannotate;

    private SagePostSnfEffApplication(final String inputVcf) {

        this.reannotate = new Reannotate(HmfGenePanelSupplier.allGenesMap37(), x -> {
//            System.out.println(x);
        });
        this.inputVcf = inputVcf;
    }

    private void run() {
        final VCFFileReader inputReader = new VCFFileReader(new File(inputVcf), false);
        for (VariantContext context : inputReader) {
            reannotate.accept(context);
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
    public void close() throws IOException {
        reannotate.close();
    }
}
