package com.hartwig.hmftools.batchvariantanalyser;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;

import javax.xml.stream.XMLStreamException;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.predicate.VariantFilter;
import com.hartwig.hmftools.common.variant.vcf.VCFFileLoader;
import com.hartwig.hmftools.common.variant.vcf.VCFSomaticFile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class BatchVariantAnalyserApplication {

    private static final Logger LOGGER = LogManager.getLogger(BatchVariantAnalyserApplication.class);

    private static final String VCF_INPUT_DIR = "vcf_input_dir";
    private static final String OUT_CSV = "out_csv";
    private static final String VARIANTS_PER_SAMPLE = "variants_per_sample";

    public static void main(final String... args)
            throws ParseException, IOException, XMLStreamException, HartwigException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        new BatchVariantAnalyserApplication(cmd.getOptionValue(VCF_INPUT_DIR), cmd.getOptionValue(OUT_CSV),
                Integer.parseInt(cmd.getOptionValue(VARIANTS_PER_SAMPLE))).run();
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();

        options.addOption(VCF_INPUT_DIR, true,
                "This path should contain a list of consensus-flagged or filtered VCFs.");
        options.addOption(OUT_CSV, true, "The list of variants will be written to this CSV.");
        options.addOption(VARIANTS_PER_SAMPLE, true, "Number of variants to select per sample");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private final String vcfInputDir;
    @NotNull
    private final String outCsv;
    private final int variantsPerSample;

    private BatchVariantAnalyserApplication(@NotNull final String vcfInputDir, @NotNull final String outCsv,
            final int variantsPerSample) {
        this.vcfInputDir = vcfInputDir;
        this.outCsv = outCsv;
        this.variantsPerSample = variantsPerSample;
    }

    private void run() throws IOException, HartwigException {
        final Multimap<String, SomaticVariant> variantsMaps = HashMultimap.create();
        final Random rand = new Random();
        for (final Path file : Files.list(new File(vcfInputDir).toPath()).collect(Collectors.toList())) {
            LOGGER.info("Processing " + file.toFile().getName());
            final VCFSomaticFile variantFile = VCFFileLoader.loadSomaticVCF(file.toFile().getPath());
            final List<SomaticVariant> variants = VariantFilter.passOnly(variantFile.variants());

            for (int i = 0; i < variantsPerSample; i++) {
                variantsMaps.put(variantFile.sample(), variants.get(rand.nextInt(variants.size())));
            }
        }
        final List<String> lines = Lists.newArrayList();
        for (final Map.Entry<String, SomaticVariant> entry : variantsMaps.entries()) {
            final SomaticVariant variant = entry.getValue();
            lines.add(
                    entry.getKey() + "," + variant.chromosome() + "," + variant.position() + "," + variant.ref() + ","
                            + variant.alt() + "," + variant.alleleReadCount() + "," + variant.totalReadCount());
        }
        Files.write(new File(outCsv).toPath(), lines);
        LOGGER.info("Written " + lines.size() + " variants to " + outCsv);
    }
}
