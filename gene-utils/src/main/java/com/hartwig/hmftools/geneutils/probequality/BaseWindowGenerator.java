package com.hartwig.hmftools.geneutils.probequality;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.deriveRefGenomeVersion;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.SpecificRegions;

// Partitions the genome into small windows of bases to be analysed individually.
public class BaseWindowGenerator
{
    private final RefGenomeSource mRefGenome;
    private final SpecificRegions mSpecificRegions;
    // Genome is partitioned into windows of this many bases.
    private final int mBaseWindowLength;
    // Windows are generated consecutively with this number of bases between the start positions.
    // < mBaseWindowLength implies overlapping windows.
    private final int mBaseWindowSpacing;
    // How many windows to yield at once.
    private final int mBatchSize;

    public BaseWindowGenerator(
            final RefGenomeSource mRefGenome, final SpecificRegions mSpecificRegions,
            final int mBaseWindowLength, final int mBaseWindowSpacing, final int mBatchSize)
    {
        this.mRefGenome = mRefGenome;
        this.mSpecificRegions = mSpecificRegions;
        if (mBaseWindowLength < 1) {
            throw new RuntimeException("mBaseWindowLength must be >= 1");
        }
        this.mBaseWindowLength = mBaseWindowLength;
        if (mBaseWindowSpacing < 1) {
            throw new RuntimeException("mBaseWindowSpacing must be >= 1");
        }
        this.mBaseWindowSpacing = mBaseWindowSpacing;
        if (mBatchSize < 1) {
            throw new RuntimeException("mBatch size must be >= 1");
        }
        this.mBatchSize = mBatchSize;
    }

    public record BaseWindow(ChrBaseRegion region, byte[] sequence) {}

    public Stream<List<BaseWindow>> createBaseWindowBatches() {
        return createBaseWindowRegionBatches()
                .map(regions ->
                        regions.stream().map(region ->
                                new BaseWindow(region, mRefGenome.getBases(region.chromosome(), region.start(), region.end()))
                        ).toList());
    }

    private Stream<List<ChrBaseRegion>> createBaseWindowRegionBatches() {
        return partitionRegionsIntoBatches(createBaseWindowRegions());
    }

    private Stream<ChrBaseRegion> createBaseWindowRegions() {
        GU_LOGGER.info("Creating base window region stream");
        RefGenomeVersion refGenomeVersion = deriveRefGenomeVersion(mRefGenome);
        RefGenomeCoordinates coordinates = refGenomeVersion.is37() ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;
        return Arrays.stream(HumanChromosome.values())
                .map(chr -> refGenomeVersion.versionedChromosome(chr.toString()))
                .filter(mSpecificRegions::includeChromosome)
                .flatMap(chr -> createBaseWindowRegions(chr, coordinates.length(chr)));
    }

    private Stream<ChrBaseRegion> createBaseWindowRegions(final String chromosome, final int chromosomeLength) {
        GU_LOGGER.debug("Creating base window stream for chromosome: {}", chromosome);
        // Generate the windows from the configured specific regions, rather than enumerating all windows and filtering, for performance.
        List<ChrBaseRegion> regions = mSpecificRegions.Regions.stream()
                .filter(region -> region.chromosome().equals(chromosome)).collect(Collectors.toList());
        if (regions.isEmpty()) {
            regions.add(new ChrBaseRegion(chromosome, 1, chromosomeLength - 1));
        }
        return regions.stream()
                .flatMap(this::createBaseWindowRegions)
                // Might be a better way to handle the end of chromosomes, for now exclude them if the window doesn't line up.
                .takeWhile(region -> region.end() <= chromosomeLength);
    }

    // Create base window regions which fully cover the specified region.
    private Stream<ChrBaseRegion> createBaseWindowRegions(final ChrBaseRegion region) {
        GU_LOGGER.debug("Creating base window stream for region: {}", region);
        // First window may start before the start of the specified region.
        int initial = baseWindowStartCoveringPosition(region.start());
        // Last window could extend past the end of the specified region.
        return IntStream.iterate(initial, start -> start <= region.end(), start -> start + mBaseWindowSpacing)
                .mapToObj(start -> new ChrBaseRegion(region.chromosome(), start, start + mBaseWindowLength - 1));
    }

    // Finds the start position of the nearest base window which covers the given position.
    private int baseWindowStartCoveringPosition(final int position) {
        int position0Idx = position - 1;
        int mod = position0Idx % mBaseWindowSpacing;
        return max(position - mod, 1);
    }

    // Batches the stream into lists of fixed size for immediate processing.
    private Stream<List<ChrBaseRegion>> partitionRegionsIntoBatches(Stream<ChrBaseRegion> regions) {
        Iterator<ChrBaseRegion> iterator = regions.iterator();
        return Stream.generate(() -> {
            List<ChrBaseRegion> batch = new ArrayList<>(mBatchSize);
            for (int i = 0; i < mBatchSize && iterator.hasNext(); i++) {
                batch.add(iterator.next());
            }
            return batch;
        }).takeWhile(b -> !b.isEmpty());
    }

    // Checks if a base sequence is "normal" for the purposes of this analysis.
    // We don't care to process regions other than regular base data.
    public static boolean isSequenceNormal(final byte[] sequence) {
        for (byte b : sequence) {
            switch (b) {
                case 'A':
                case 'C':
                case 'G':
                case 'T':
                case 'a':
                case 'c':
                case 'g':
                case 't':
                    continue;
                default:
                    return false;
            }
        }
        return true;
    }
}
