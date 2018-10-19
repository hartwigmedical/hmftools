package com.hartwig.hmftools.purple.config;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.ChromosomeLengthFile;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositionFactory;
import com.hartwig.hmftools.common.refgenome.RefGenome;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface RefGenomeConfig {

    Logger LOGGER = LogManager.getLogger(RefGenomeConfig.class);

    String REF_GENOME = "ref_genome";

    static void addOptions(@NotNull Options options) {
        options.addOption(REF_GENOME,
                true,
                "Reference genome to use. Will attempt to detect using cobalt chromosome lengths, otherwise must be either \"hg19\" or \"hg38\".");
    }

    @NotNull
    Map<Chromosome, GenomePosition> length();

    @NotNull
    Map<Chromosome, GenomePosition> centromere();

    @NotNull
    static RefGenomeConfig createRefGenomeConfig(@NotNull CommandLine cmd, @NotNull final String tumorSample,
            @NotNull final String cobaltDirectory) throws IOException, ParseException {
        final String chrLengthFile = ChromosomeLengthFile.generateFilename(cobaltDirectory, tumorSample);
        if (!new File(chrLengthFile).exists()) {
            throw new ParseException("Unable to locate cobalt chromosome length file " + chrLengthFile);
        }

        final Map<Chromosome, GenomePosition> lengthPositions = ChromosomeLengthFile.read(chrLengthFile)
                .stream()
                .filter(x -> HumanChromosome.contains(x.chromosome()))
                .collect(Collectors.toMap(x -> HumanChromosome.fromString(x.chromosome()),
                        item -> GenomePositionFactory.create(item.chromosome(), item.length())));

        final Map<Chromosome, Long> lengths = lengthPositions.values()
                .stream()
                .collect(Collectors.toMap(x -> HumanChromosome.fromString(x.chromosome()), GenomePosition::position));

        final Optional<RefGenome> automaticallyDetectedRefGenome = RefGenome.fromLengths(lengths);

        final RefGenome refGenome;
        if (cmd.hasOption(REF_GENOME)) {
            try {
                refGenome = RefGenome.valueOf(cmd.getOptionValue(REF_GENOME).toUpperCase());
                LOGGER.info("Using parameter supplied ref genome: {}", refGenome);
            } catch (Exception e) {
                throw new ParseException(
                        "Unknown ref genome " + cmd.getOptionValue(REF_GENOME) + ". Must be either \"hg19\" or \"hg38\".");
            }
        } else if (automaticallyDetectedRefGenome.isPresent()) {
            refGenome = automaticallyDetectedRefGenome.get();
            LOGGER.info("Detected ref genome: {}", refGenome);
        } else {
            throw new ParseException("Unable to detect ref genome. Please specify " + cmd.getOptionValue(REF_GENOME)
                    + " parameter as one of \"hg19\" or \"hg38\". ");
        }

        return ImmutableRefGenomeConfig.builder().length(lengthPositions).centromere(centromeres(refGenome, lengthPositions)).build();

    }

    @NotNull
    static Map<Chromosome, GenomePosition> centromeres(@NotNull final RefGenome refGenome,
            @NotNull final Map<Chromosome, GenomePosition> lengths) {
        final Map<Chromosome, Long> centromeres = refGenome.centromeres();
        final Map<Chromosome, GenomePosition> result = Maps.newHashMap();

        for (Map.Entry<Chromosome, GenomePosition> entry : lengths.entrySet()) {
            final Chromosome chromosome = entry.getKey();
            final String contig = entry.getValue().chromosome();
            if (centromeres.containsKey(chromosome)) {
                result.put(chromosome, GenomePositionFactory.create(contig, centromeres.get(chromosome)));
            }

        }

        return result;
    }
}
