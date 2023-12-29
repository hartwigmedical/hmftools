package com.hartwig.hmftools.esvee.read;

import static java.lang.Math.max;

import static com.hartwig.hmftools.esvee.SvConstants.BAM_READ_JUNCTION_BUFFER;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.common.RegionOfInterest;
import com.hartwig.hmftools.esvee.util.MultiMapCollector;

import org.apache.commons.lang3.tuple.Pair;

public interface SAMSource extends AutoCloseable
{
    Stream<Read> streamReadsContaining(final String chromosome, final int startPosition, final int endPosition);

    default List<Read> findReadsContaining(final ChrBaseRegion region)
    {
        return streamReadsContaining(region.chromosome(), region.start(), region.end()).collect(Collectors.toList());
    }

    default List<Read> findReadsContaining(final String chromosome, final int startPosition, final int endPosition)
    {
        return streamReadsContaining(chromosome, startPosition, endPosition).collect(Collectors.toList());
    }

    default List<Read> findReadsNear(final String chromosome, final int position)
    {
        int bamReadStart = max(position - BAM_READ_JUNCTION_BUFFER, 1);
        int bamReadEnd = position + BAM_READ_JUNCTION_BUFFER; // CHECK: could go over end of chromosome and cause an error
        return streamReadsContaining(chromosome, bamReadStart, bamReadEnd).collect(Collectors.toList());
    }

    // returns tuples showing the original record on the left, and the mate on the right.
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
                .map(mateRegion -> streamReadsContaining(mateRegion.Chromosome, mateRegion.start(), mateRegion.end())
                        .flatMap(record -> recordsByName.getOrDefault(record.getName(), List.of()).stream()
                                .filter(match -> match.isFirstOfPair() != record.isFirstOfPair())
                                .map(match -> Pair.of(match, record))))
                .reduce(Stream::concat).orElse(Stream.of());
    }

    @Override
    void close();
}
