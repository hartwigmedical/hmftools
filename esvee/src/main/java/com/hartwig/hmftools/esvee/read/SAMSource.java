package com.hartwig.hmftools.esvee.read;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.esvee.common.RegionOfInterest;
import com.hartwig.hmftools.esvee.util.MultiMapCollector;

import org.apache.commons.lang3.tuple.Pair;

public interface SAMSource extends AutoCloseable
{
    Stream<Read> streamReadsContaining(final String chromosome, final int startPosition, final int endPosition);

    default List<Read> findReadsContaining(final RegionOfInterest region)
    {
        return findReadsContaining(region.Chromosome, region.Start, region.End);
    }

    default List<Read> findReadsContaining(final String chromosome, final int startPosition, final int endPosition)
    {
        return streamReadsContaining(chromosome, startPosition, endPosition).collect(Collectors.toList());
    }

    default Stream<Read> streamReadsNear(final String chromosome, final int position)
    {
        return streamReadsContaining(chromosome, Math.max(position - 1000, 1), position + 1000);
    }

    default List<Read> findReadsNear(final String chromosome, final int position)
    {
        return streamReadsNear(chromosome, position).collect(Collectors.toList());
    }

    /**
     * @return tuples showing the original record on the left, and the mate on the right.
     */
    default Stream<Pair<Read, Read>> streamMates(final Collection<Read> reads)
    {
        final Stream<RegionOfInterest> mateRegionsUnmerged = reads.stream()
                .map(record -> record.isMateMapped()
                        ? new RegionOfInterest(record.getMateChromosome(), record.getMateAlignmentStart(),
                        record.getMateAlignmentStart() + record.getLength())
                        : new RegionOfInterest(record.getChromosome(), record.getAlignmentStart(),
                                record.getAlignmentStart() + record.getLength()));
        final List<RegionOfInterest> mappedMateRegions = RegionOfInterest.tryMerge(mateRegionsUnmerged::iterator);

        final Map<String, List<Read>> recordsByName = reads.stream()
                .collect(MultiMapCollector.keyed(Read::getName));

        return mappedMateRegions.stream()
                .map(mateRegion -> streamReadsContaining(mateRegion.Chromosome, mateRegion.Start, mateRegion.End)
                        .flatMap(record -> recordsByName.getOrDefault(record.getName(), List.of()).stream()
                                .filter(match -> match.isFirstOfPair() != record.isFirstOfPair())
                                .map(match -> Pair.of(match, record))))
                .reduce(Stream::concat).orElse(Stream.of());
    }

    default List<Pair<Read, Read>> findMates(final Collection<Read> reads)
    {
        return streamMates(reads).collect(Collectors.toList());
    }

    @Override // No exceptions
    void close();
}
