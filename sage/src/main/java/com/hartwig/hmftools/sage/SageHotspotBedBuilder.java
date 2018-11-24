package com.hartwig.hmftools.sage;

import static java.util.Collections.emptyList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.collect.Multimaps;
import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodSupplier;
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.hotspot.VariantHotspotFile;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionBuilder;
import com.hartwig.hmftools.common.region.HmfExonRegion;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class SageHotspotBedBuilder {

    private static final Logger LOGGER = LogManager.getLogger(SageHotspotBedBuilder.class);

    private static final String OUT_PATH = "out";
    private static final String HOTSPOT = "hotspot";
    private static final String SEPARATOR = "\t";

    public static void main(String[] args) throws IOException, ParseException {

        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String outputFilePath = cmd.getOptionValue(OUT_PATH);
        final String hotspotPath = cmd.getOptionValue(HOTSPOT);

        if (outputFilePath == null || hotspotPath == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("SageHotspotBedBuilder", options);
            System.exit(1);
        }

        final Set<String> driverGenes = DndsDriverGeneLikelihoodSupplier.oncoLikelihood().keySet();
        LOGGER.info("Loaded {} oncogenes", driverGenes.size());

        final ListMultimap<Chromosome, VariantHotspot> hotspots = VariantHotspotFile.read(hotspotPath);
        LOGGER.info("Loaded {} hotspots from {}", hotspots.values().size(), hotspotPath);

        final ListMultimap<Chromosome, HmfTranscriptRegion> geneRegions = Multimaps.fromRegions(HmfGenePanelSupplier.allGeneList37()
                .stream()
                .filter(x -> driverGenes.contains(x.gene()))
                .collect(Collectors.toList()));

        LOGGER.info("Merging oncogene coding regions with known hotspot locations");
        final List<String> bedResult = Lists.newArrayList();
        for (HumanChromosome chromosome : HumanChromosome.values()) {
            List<GenomeRegion> result = Lists.newArrayList();
            final List<HmfTranscriptRegion> chromosomeRegions =
                    geneRegions.containsKey(chromosome) ? geneRegions.get(chromosome) : emptyList();
            for (HmfTranscriptRegion gene : chromosomeRegions) {
                result.addAll(addCodingRegions(gene));
            }

            final List<VariantHotspot> chromosomeHotspots = hotspots.containsKey(chromosome) ? hotspots.get(chromosome) : emptyList();
            for (VariantHotspot hotspot : chromosomeHotspots) {
                result.addAll(addVariantHotspot(hotspot));
            }

            result.stream().map(SageHotspotBedBuilder::toBedFormat).forEach(bedResult::add);
        }

        LOGGER.info("Writing {} regions to bed file {}", bedResult.size(), outputFilePath);
        Files.write(new File(outputFilePath).toPath(), bedResult);
    }

    @NotNull
    private static List<GenomeRegion> addCodingRegions(@NotNull final HmfTranscriptRegion gene) {
        final GenomeRegionBuilder builder = new GenomeRegionBuilder(gene.chromosome());
        for (HmfExonRegion exon : gene.exome()) {
            for (long position = exon.start(); position <= exon.end(); position++) {
                if (position >= gene.codingStart() && position <= gene.codingEnd()) {
                    builder.addPosition(position);
                }
            }
        }

        return builder.build();
    }

    @NotNull
    static List<GenomeRegion> addVariantHotspot(@NotNull final VariantHotspot hotspot) {
        final GenomeRegionBuilder builder = new GenomeRegionBuilder(hotspot.chromosome());
        for (int i = 0; i < Math.min(hotspot.ref().length(), hotspot.alt().length()); i++) {
            builder.addPosition(hotspot.position() + i);
        }

        return builder.build();
    }

    @NotNull
    private static String toBedFormat(@NotNull final GenomeRegion exon) {
        return exon.chromosome() + SEPARATOR + (exon.start() - 1) + SEPARATOR + exon.end();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(OUT_PATH, true, "Bed output file.");
        options.addOption(HOTSPOT, true, "Hotspot input file.");
        return options;
    }
}
