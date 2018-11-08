package com.hartwig.hmftools.hotspotcaller;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.hotspot.HotspotEvidence;
import com.hartwig.hmftools.common.hotspot.HotspotEvidenceFactory;
import com.hartwig.hmftools.common.hotspot.HotspotEvidenceVCF;
import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.hotspot.VariantHotspotFile;
import com.hartwig.hmftools.common.pileup.Pileup;
import com.hartwig.hmftools.common.pileup.PileupFile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class HotspotCallerApplication {

    private static final Logger LOGGER = LogManager.getLogger(HotspotCallerApplication.class);

    private static final String OUT_PATH = "out";
    private static final String HOTSPOT = "hotspot";
    private static final String REFERENCE = "reference";
    private static final String REFERENCE_PILEUP = "reference_pileup";
    private static final String TUMOR = "tumor";
    private static final String TUMOR_PILEUP = "tumor_pileup";

    public static void main(String[] args) throws IOException, ParseException {

        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String outputFilePath = cmd.getOptionValue(OUT_PATH);
        final String hotspotPath = cmd.getOptionValue(HOTSPOT);
        final String tumorSample = cmd.getOptionValue(TUMOR);
        final String tumorPileupPath = cmd.getOptionValue(TUMOR_PILEUP);
        final String referenceSample = cmd.getOptionValue(REFERENCE);
        final String referencePileupPath = cmd.getOptionValue(REFERENCE_PILEUP);

        if (outputFilePath == null || hotspotPath == null || tumorSample == null || tumorPileupPath == null || referenceSample == null
                || referencePileupPath == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("HotspotCallerApplication", options);
            System.exit(1);
        }

        LOGGER.info("Loading hotspots from {}", hotspotPath);
        final ListMultimap<Chromosome, VariantHotspot> hotspots = VariantHotspotFile.read(hotspotPath);

        LOGGER.info("Loading reference pileup from {}", referencePileupPath);
        final List<Pileup> normalPileup = PileupFile.read(referencePileupPath);

        LOGGER.info("Loading tumor pileup from {}", tumorPileupPath);
        final List<Pileup> tumorPileup = PileupFile.read(tumorPileupPath);

        LOGGER.info("Examining evidence.");
        final HotspotEvidenceFactory evidenceFactory = new HotspotEvidenceFactory(hotspots);
        final List<HotspotEvidence> evidence = evidenceFactory.evidence(tumorPileup, normalPileup);

        LOGGER.info("Writing output to {}", outputFilePath);
        new HotspotEvidenceVCF(referenceSample, tumorSample).write(outputFilePath, evidence);

        LOGGER.info("Complete");
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(OUT_PATH, true, "Output file.");
        options.addOption(HOTSPOT, true, "Hotspot input file.");
        options.addOption(TUMOR, true, "Tumor sample name.");
        options.addOption(TUMOR_PILEUP, true, "Tumor pileup path.");
        options.addOption(REFERENCE, true, "Reference sample name.");
        options.addOption(REFERENCE_PILEUP, true, "Reference pileup oath.");
        return options;
    }
}
