package com.hartwig.hmftools.purple.config;

import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.common.refgenome.RefGenome;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.region.ImmutableHmfExonRegion;
import com.hartwig.hmftools.common.region.ImmutableHmfTranscriptRegion;

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
public interface RefGenomeData {

    Logger LOGGER = LogManager.getLogger(RefGenomeData.class);

    String REF_GENOME = "ref_genome";

    static void addOptions(@NotNull Options options) {
        options.addOption(REF_GENOME,
                true,
                "Reference genome to use. Will attempt to detect using cobalt chromosome lengths, otherwise must be either \"hg19\" or \"hg38\".");
    }

    @NotNull
    RefGenome refRegome();

    @NotNull
    Map<Chromosome, GenomePosition> length();

    @NotNull
    Map<Chromosome, GenomePosition> centromere();

    @NotNull
    List<HmfTranscriptRegion> genePanel();

    @NotNull
    static RefGenomeData createRefGenomeConfig(@NotNull CommandLine cmd, @NotNull CobaltData cobaltData) throws ParseException {

        if (cobaltData.chromosomeLengthsEstimated() && !cmd.hasOption(REF_GENOME)) {
            throw new ParseException(
                    "Cobalt chromosome information unavailable. You must specify -" + REF_GENOME + " as either \"hg19\" or \"hg38\".");
        }

        final Map<Chromosome, GenomePosition> lengthPositions = cobaltData.chromosomeLengths();
        final Map<Chromosome, Long> lengths = asLongs(lengthPositions);
        final Optional<RefGenome> automaticallyDetectedRefGenome = RefGenome.fromLengths(lengths);
        final Map<Chromosome, String> contigMap =
                lengthPositions.entrySet().stream().collect(Collectors.toMap(Map.Entry::getKey, x -> x.getValue().chromosome()));

        final Optional<RefGenome> suppliedGenome;
        if (cmd.hasOption(REF_GENOME)) {
            try {
                suppliedGenome = Optional.of(RefGenome.valueOf(cmd.getOptionValue(REF_GENOME).toUpperCase()));
            } catch (Exception e) {
                throw new ParseException("Unknown ref genome " + cmd.getOptionValue(REF_GENOME) + ". Must be either \"hg19\" or \"hg38\".");
            }
        } else {
            suppliedGenome = Optional.empty();
        }

        final RefGenome refGenome;
        if (automaticallyDetectedRefGenome.isPresent()) {
            refGenome = automaticallyDetectedRefGenome.get();
            LOGGER.info("Detected ref genome: {}", refGenome);
            if (suppliedGenome.isPresent() && !suppliedGenome.get().equals(automaticallyDetectedRefGenome.get())) {
                throw new ParseException("Parameter " + REF_GENOME + " " + suppliedGenome.get() + " does not match detected ref genome.");
            }
        } else if (suppliedGenome.isPresent()) {
            refGenome = suppliedGenome.get();
            LOGGER.info("Using {} parameter: {}", REF_GENOME, refGenome);
        } else {
            throw new ParseException("Unable to detect ref genome. Please specify " + cmd.getOptionValue(REF_GENOME)
                    + " parameter as one of \"hg19\" or \"hg38\". ");
        }

        final List<HmfTranscriptRegion> rawGenePanel =
                refGenome == RefGenome.HG38 ? HmfGenePanelSupplier.allGeneList38() : HmfGenePanelSupplier.allGeneList37();

        final Function<String, String> chromosomeToContig = s -> contigMap.get(HumanChromosome.fromString(s));

        final List<HmfTranscriptRegion> genePanel = rawGenePanel.stream()
                .filter(x -> HumanChromosome.contains(x.chromosome()) && contigMap.containsKey(HumanChromosome.fromString(x.chromosome())))
                .map(x -> updateContig(x, chromosomeToContig))
                .collect(Collectors.toList());

        return ImmutableRefGenomeData.builder()
                .length(toPosition(refGenome.lengths(), contigMap))
                .centromere(toPosition(refGenome.centromeres(), contigMap))
                .refRegome(refGenome)
                .genePanel(genePanel)
                .build();
    }

    @NotNull
    static HmfTranscriptRegion updateContig(HmfTranscriptRegion region, Function<String, String> fixContigFunction) {
        final String correctContig = fixContigFunction.apply(region.chromosome());

        return ImmutableHmfTranscriptRegion.builder()
                .from(region)
                .chromosome(correctContig)
                .exome(region.exome()
                        .stream()
                        .map(x -> ImmutableHmfExonRegion.builder().from(x).chromosome(correctContig).build())
                        .collect(Collectors.toList()))
                .build();
    }

    @NotNull
    static Map<Chromosome, GenomePosition> toPosition(@NotNull final Map<Chromosome, Long> longs,
            @NotNull final Map<Chromosome, String> contigMap) {
        final Map<Chromosome, GenomePosition> result = Maps.newHashMap();

        for (Map.Entry<Chromosome, String> entry : contigMap.entrySet()) {
            final Chromosome chromosome = entry.getKey();
            final String contig = entry.getValue();
            if (longs.containsKey(chromosome)) {
                result.put(chromosome, GenomePositions.create(contig, longs.get(chromosome)));
            }

        }

        return result;
    }

    static Map<Chromosome, Long> asLongs(@NotNull final Map<Chromosome, GenomePosition> positions) {
        return positions.values()
                .stream()
                .collect(Collectors.toMap(x -> HumanChromosome.fromString(x.chromosome()), GenomePosition::position));
    }

}
