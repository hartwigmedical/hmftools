package com.hartwig.hmftools.geneutils.probequality;

import static java.lang.Math.max;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.SpecificRegions;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

// Partitions the genome into small windows of bases to be analysed individually.
public class BaseWindowGenerator
{
    private final RefGenomeInterface mRefGenome;
    private final SpecificRegions mSpecificRegions;
    // Genome is partitioned into windows of this many bases.
    private final int mBaseWindowLength;
    // Windows are generated consecutively with this number of bases between the start positions.
    // < mBaseWindowLength implies overlapping windows.
    private final int mBaseWindowSpacing;
    // How many windows to yield at once.
    private final int mBatchSize;

    private static final Logger LOGGER = LogManager.getLogger(BaseWindowGenerator.class);

    public BaseWindowGenerator(
            final RefGenomeInterface mRefGenome,
            final SpecificRegions mSpecificRegions,
            final int mBaseWindowLength, final int mBaseWindowSpacing, final int mBatchSize)
    {
        this.mRefGenome = mRefGenome;
        this.mSpecificRegions = mSpecificRegions;
        if(mBaseWindowLength < 1)
        {
            throw new IllegalArgumentException("mBaseWindowLength must be >= 1");
        }
        this.mBaseWindowLength = mBaseWindowLength;
        if(mBaseWindowSpacing < 1)
        {
            throw new IllegalArgumentException("mBaseWindowSpacing must be >= 1");
        }
        this.mBaseWindowSpacing = mBaseWindowSpacing;
        if(mBatchSize < 1)
        {
            throw new IllegalArgumentException("mBatch size must be >= 1");
        }
        this.mBatchSize = mBatchSize;
    }

    public record BaseWindow(ChrBaseRegion region, byte[] sequence)
    {
        public BaseWindow
        {
            if(region.baseLength() != sequence.length)
            {
                throw new IllegalArgumentException("BaseWindow region and sequence must be same length");
            }
        }

        @Override
        public boolean equals(final Object obj)
        {
            if(obj instanceof BaseWindow)
            {
                return region.equals(((BaseWindow) obj).region) && Arrays.equals(sequence, ((BaseWindow) obj).sequence);
            }
            else
            {
                return false;
            }
        }

        @NotNull
        @Override
        public String toString()
        {
            return String.format("%s(%s)", region.toString(), new String(sequence));
        }
    }

    public Stream<List<BaseWindow>> createBaseWindowBatches()
    {
        return createBaseWindowRegionBatches()
                .map(regions ->
                        regions.stream().map(region ->
                                new BaseWindow(region, mRefGenome.getBases(region.chromosome(), region.start(), region.end()))
                        ).toList());
    }

    public Stream<List<ChrBaseRegion>> createBaseWindowRegionBatches()
    {
        return batchRegions(createBaseWindowRegions());
    }

    public Stream<ChrBaseRegion> createBaseWindowRegions()
    {
        LOGGER.debug("Creating base window region stream");

        return mRefGenome.chromosomeLengths().keySet().stream().sorted().flatMap(this::createBaseWindowRegions);
    }

    public Stream<ChrBaseRegion> createBaseWindowRegions(final String chromosome)
    {
        LOGGER.debug("Creating base window stream for chromosome: {}", chromosome);

        // Note if any chromosomes are specified then only those chromosome are examined, regardless of other filters.
        if(!mSpecificRegions.includeChromosome(chromosome))
        {
            return Stream.empty();
        }

        int chromosomeLength = mRefGenome.getChromosomeLength(chromosome);

        // Generate the windows from the configured specific regions, rather than enumerating all windows and filtering, for performance.
        ArrayList<ChrBaseRegion> regions = mSpecificRegions.Regions.stream()
                .filter(region -> region.chromosome().equals(chromosome))
                .collect(Collectors.toCollection(ArrayList::new));
        // If no specific regions then add the whole chromosome.
        if(regions.isEmpty())
        {
            regions.add(new ChrBaseRegion(chromosome, 1, chromosomeLength - 1));
        }

        // Ensure the regions are within the chromosome bounds.
        regions.forEach(region ->
        {
            if(!(region.start() >= 1 && region.end() <= chromosomeLength))
            {
                throw new IllegalArgumentException("Specified region is out of bounds");
            }
        });

        // Merge overlapping regions otherwise we could generate some windows more than once.
        ChrBaseRegion.checkMergeOverlaps(regions, true);

        return regions.stream()
                .flatMap(this::createBaseWindowRegions)
                // Might be a better way to handle the end of chromosomes, for now exclude them if the window doesn't line up.
                // Note this could exclude part of a configured specific region, which may be surprising to the user.
                .takeWhile(region -> region.end() <= chromosomeLength);
    }

    // Create base window regions which fully cover the specified region.
    public Stream<ChrBaseRegion> createBaseWindowRegions(final ChrBaseRegion region)
    {
        LOGGER.debug("Creating base window stream for region: {}", region);
        if(region.baseLength() <= 1)
        {
            return Stream.empty();
        }
        // First window may start before the start of the specified region.
        int initial = baseWindowStartCoveringPosition(region.start());
        // Last window could extend past the end of the specified region.
        return IntStream.iterate(initial, start -> start <= region.end(), start -> start + mBaseWindowSpacing)
                .mapToObj(start -> new ChrBaseRegion(region.chromosome(), start, start + mBaseWindowLength - 1));
    }

    // Finds the start position of the nearest base window which covers the given position.
    public int baseWindowStartCoveringPosition(final int position)
    {
        if(position < 1)
        {
            throw new IllegalArgumentException("position must be 1-indexed");
        }
        int position0Idx = position - 1;
        int mod = position0Idx % mBaseWindowSpacing;
        return max(position - mod, 1);
    }

    // Batches the stream into lists of fixed size for immediate processing.
    public Stream<List<ChrBaseRegion>> batchRegions(Stream<ChrBaseRegion> regions)
    {
        Iterator<ChrBaseRegion> iterator = regions.iterator();
        return Stream.generate(() ->
        {
            List<ChrBaseRegion> batch = new ArrayList<>(mBatchSize);
            for(int i = 0; i < mBatchSize && iterator.hasNext(); i++)
            {
                batch.add(iterator.next());
            }
            return batch;
        }).takeWhile(b -> !b.isEmpty());
    }

    // Checks if a base sequence is "normal" for the purposes of this analysis.
    // We don't care to process regions other than regular base data.
    public static boolean isSequenceNormal(final byte[] sequence)
    {
        for(byte b : sequence)
        {
            switch(b)
            {
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
