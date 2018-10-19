package com.hartwig.hmftools.purple.segment;

import static com.hartwig.hmftools.common.position.GenomePositions.union;

import java.io.IOException;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.centromeres.Centromeres;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.pcf.ImmutablePCFPosition;
import com.hartwig.hmftools.common.pcf.PCFFile;
import com.hartwig.hmftools.common.pcf.PCFPosition;
import com.hartwig.hmftools.common.pcf.PCFSource;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.window.Window;
import com.hartwig.hmftools.purple.config.CommonConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class PCFPositionsSupplier {
    private static final Logger LOGGER = LogManager.getLogger(PCFPositionsSupplier.class);

    @NotNull
    public static Multimap<String, PCFPosition> createPositions(@NotNull final CommonConfig config) throws IOException {
        final String referenceFile = PCFFile.generateRatioFilename(config.cobaltDirectory(), config.refSample());
        LOGGER.info("Loading reference ratio PCF segments from {}", referenceFile);
        final Multimap<String, PCFPosition> referenceBreakPoint =
                PCFFile.readPositions(config.windowSize(), PCFSource.REFERENCE_RATIO, referenceFile);

        final String tumorFile = PCFFile.generateRatioFilename(config.cobaltDirectory(), config.tumorSample());
        LOGGER.info("Loading tumor ratio PCF segments from {}", tumorFile);
        final Multimap<String, PCFPosition> tumorBreakPoints = PCFFile.readPositions(config.windowSize(), PCFSource.TUMOR_RATIO, tumorFile);

        final String amberFile = PCFFile.generateBAFFilename(config.amberDirectory(), config.tumorSample());
        LOGGER.info("Loading tumor baf PCF segments from {}", amberFile);
        final Multimap<String, PCFPosition> tumorBAF = PCFFile.readPositions(config.windowSize(), PCFSource.TUMOR_BAF, amberFile);

        final Set<String> contigs = tumorBAF.values().stream().map(GenomePosition::chromosome).collect(Collectors.toSet());

        LOGGER.info("Merging reference and tumor ratio break points");
        return union(union(union(tumorBreakPoints, referenceBreakPoint), tumorBAF), centromeres(config.windowSize(), contigs));
    }

    @NotNull
    private static Multimap<String, PCFPosition> centromeres(int windowSize, @NotNull Set<String> contigs) {
        final Multimap<String, PCFPosition> result = ArrayListMultimap.create();
        final Window window = new Window(windowSize);
        final Map<Chromosome, Long> centromeres = Centromeres.hg19();
        for (String contig : contigs) {
            final Chromosome chromosome = HumanChromosome.fromString(contig);
            final Long centromere = centromeres.get(chromosome);
            if (centromere != null) {
                result.put(contig, create(contig, window.start(centromere)));
            }

        }

        return result;
    }

    @NotNull
    private static PCFPosition create(@NotNull final String chromosome, long position) {
        return ImmutablePCFPosition.builder()
                .chromosome(chromosome)
                .position(position)
                .source(PCFSource.REFERENCE_RATIO)
                .minPosition(position)
                .maxPosition(position)
                .build();
    }
}
