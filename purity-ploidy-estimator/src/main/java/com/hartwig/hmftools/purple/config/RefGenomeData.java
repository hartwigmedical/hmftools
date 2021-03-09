package com.hartwig.hmftools.purple.config;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.ChromosomeLength;
import com.hartwig.hmftools.common.genome.chromosome.ChromosomeLengthFactory;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface RefGenomeData {

    Logger LOGGER = LogManager.getLogger(RefGenomeData.class);

    String REF_GENOME = "ref_genome";

    static void addOptions(@NotNull Options options) {
        options.addOption(REF_GENOME, true, "Path to the ref genome fasta file.");
    }

    boolean isHg38();

    @NotNull
    String refGenome();

    @NotNull
    Map<Chromosome, GenomePosition> length();

    @NotNull
    Map<Chromosome, GenomePosition> centromere();

    @NotNull
    List<HmfTranscriptRegion> hmfTranscripts();

    @NotNull
    static RefGenomeData createRefGenomeConfig(@NotNull CommandLine cmd) throws ParseException, IOException {

        if (!cmd.hasOption(REF_GENOME)) {
            throw new ParseException(REF_GENOME + " is a mandatory argument");
        }

        final String refGenomePath = cmd.getOptionValue(REF_GENOME);
        final Map<Chromosome, GenomePosition> lengthPositions;
        try (final IndexedFastaSequenceFile indexedFastaSequenceFile = new IndexedFastaSequenceFile(new File(refGenomePath))) {
            SAMSequenceDictionary sequenceDictionary = indexedFastaSequenceFile.getSequenceDictionary();
            if (sequenceDictionary == null) {
                throw new ParseException("Supplied ref genome must have associated sequence dictionary");
            }

            lengthPositions = fromLengths(ChromosomeLengthFactory.create(indexedFastaSequenceFile.getSequenceDictionary()));
        }

        final GenomePosition chr1Length = lengthPositions.get(HumanChromosome._1);
        final RefGenomeCoordinates refGenome;
        if (chr1Length != null && chr1Length.position() == RefGenomeCoordinates.COORDS_38.lengths().get(HumanChromosome._1)) {
            refGenome = RefGenomeCoordinates.COORDS_38;
        } else {
            refGenome = RefGenomeCoordinates.COORDS_37;
        }
        LOGGER.info("Using ref genome: {}", refGenome);

        final Map<Chromosome, String> contigMap =
                lengthPositions.entrySet().stream().collect(Collectors.toMap(Map.Entry::getKey, x -> x.getValue().chromosome()));

        final List<HmfTranscriptRegion> genePanel =
                refGenome == RefGenomeCoordinates.COORDS_38 ? HmfGenePanelSupplier.allGeneList38() : HmfGenePanelSupplier.allGeneList37();

        return ImmutableRefGenomeData.builder()
                .length(toPosition(refGenome.lengths(), contigMap))
                .centromere(toPosition(refGenome.centromeres(), contigMap))
                .refGenome(refGenomePath)
                .isHg38(refGenome.equals(RefGenomeCoordinates.COORDS_38))
                .hmfTranscripts(genePanel)
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

    @NotNull
    static Map<Chromosome, GenomePosition> fromLengths(@NotNull final Collection<ChromosomeLength> lengths) {
        return lengths.stream()
                .filter(x -> HumanChromosome.contains(x.chromosome()))
                .collect(Collectors.toMap(x -> HumanChromosome.fromString(x.chromosome()),
                        item -> GenomePositions.create(item.chromosome(), item.length())));
    }

}
