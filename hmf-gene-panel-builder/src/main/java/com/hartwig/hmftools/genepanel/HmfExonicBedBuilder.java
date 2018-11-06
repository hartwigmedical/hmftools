package com.hartwig.hmftools.genepanel;

import static java.util.Collections.emptyList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.collect.Multimaps;
import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodSupplier;
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionFactory;
import com.hartwig.hmftools.common.region.HmfExonRegion;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.hotspot.VariantHotspotFile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class HmfExonicBedBuilder {

    private static final Logger LOGGER = LogManager.getLogger(HmfExonicBedBuilder.class);

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
            formatter.printHelp("HmfExonicBedBuilder", options);
            System.exit(1);
        }

        final Set<String> driverGenes = Sets.newHashSet();
        driverGenes.addAll(DndsDriverGeneLikelihoodSupplier.tsgLikelihood().keySet());
        driverGenes.addAll(DndsDriverGeneLikelihoodSupplier.oncoLikelihood().keySet());

        LOGGER.info("Loading hotspots from {}", hotspotPath);
        final ListMultimap<Chromosome, VariantHotspot> hotspots = VariantHotspotFile.read(hotspotPath);

        LOGGER.info("Loading gene definitions");
        final ListMultimap<Chromosome, HmfTranscriptRegion> geneRegions = Multimaps.fromRegions(HmfGenePanelSupplier.allGeneList37()
                .stream()
                .filter(x -> driverGenes.contains(x.gene()))
                .collect(Collectors.toList()));

        LOGGER.info("Processing");
        final List<String> bedResult = Lists.newArrayList();
        for (HumanChromosome chromosome : HumanChromosome.values()) {
            List<GenomeRegion> result = Lists.newArrayList();
            final List<HmfTranscriptRegion> chromosomeRegions =
                    geneRegions.containsKey(chromosome) ? geneRegions.get(chromosome) : emptyList();
            for (HmfTranscriptRegion gene : chromosomeRegions) {
                result = addGene(gene, result);
            }

            final List<VariantHotspot> chromosomeHotspots = hotspots.containsKey(chromosome) ? hotspots.get(chromosome) : emptyList();
            for (VariantHotspot hotspot : chromosomeHotspots) {
                result = addVariantHotspot(hotspot, result);
            }

            result.stream().map(HmfExonicBedBuilder::toBedFormat).forEach(bedResult::add);
        }

        LOGGER.info("Writing bed output {}", outputFilePath);
        Files.write(new File(outputFilePath).toPath(), bedResult);
    }

    @NotNull
    private static List<GenomeRegion> addGene(@NotNull final HmfTranscriptRegion gene, @NotNull final List<GenomeRegion> region) {
        List<GenomeRegion> result = region;
        for (HmfExonRegion exon : gene.exome()) {
            for (long position = exon.start() - 2; position <= exon.end() + 2; position++) {
                result = addPosition(exon.chromosome(), position, result);
            }
        }

        return result;
    }

    @NotNull
    static List<GenomeRegion> addVariantHotspot(@NotNull final VariantHotspot hotspot, @NotNull final List<GenomeRegion> region) {
        List<GenomeRegion> result = region;
        for (int i = 0; i < Math.min(hotspot.ref().length(), hotspot.alt().length()); i++) {
            result = addPosition(hotspot.chromosome(), hotspot.position() + i, result);
        }

        return result;
    }

    @NotNull
    static List<GenomeRegion> addPosition(@NotNull final String chromosome, final long position,
            @NotNull final List<GenomeRegion> regions) {
        GenomeRegion prev = null;

        for (int i = 0; i < regions.size(); i++) {
            GenomeRegion current = regions.get(i);
            if (position >= current.start() && position <= current.end()) {
                return regions;
            }

            if (position < current.start()) {
                // Attach to previous
                if (prev != null && position == prev.end() + 1) {
                    if (current.start() == prev.end() + 2) {
                        prev = GenomeRegionFactory.create(prev.chromosome(), prev.start(), current.end());
                        regions.set(i - 1, prev);
                        regions.remove(i);
                        return regions;
                    } else {
                        prev = GenomeRegionFactory.create(prev.chromosome(), prev.start(), prev.end() + 1);
                        regions.set(i - 1, prev);
                        return regions;
                    }
                }

                // Attach to current
                if (position == current.start() - 1) {
                    current = GenomeRegionFactory.create(current.chromosome(), current.start() - 1, current.end());
                    regions.set(i, current);
                    return regions;
                }

                // Attach between
                regions.add(i, GenomeRegionFactory.create(chromosome, position, position));
                return regions;
            }

            prev = current;
        }

        if (prev != null && position == prev.end() + 1) {
            prev = GenomeRegionFactory.create(prev.chromosome(), prev.start(), prev.end() + 1);
            regions.set(regions.size() - 1, prev);
        } else {
            regions.add(GenomeRegionFactory.create(chromosome, position, position));
        }

        return regions;
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
