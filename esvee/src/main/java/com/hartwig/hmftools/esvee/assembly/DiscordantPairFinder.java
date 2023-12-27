package com.hartwig.hmftools.esvee.assembly;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.esvee.RegionOfInterest;
import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.models.Record;
import com.hartwig.hmftools.esvee.sam.SAMSource;
import com.hartwig.hmftools.esvee.util.SizedIterable;

import org.apache.commons.lang3.tuple.Pair;

public class DiscordantPairFinder
{
    private final SAMSource mSource;

    private final int mSearchDistance;
    private final int mMinQuality;
    private final int mMinFragmentLength;

    public DiscordantPairFinder(final SAMSource source)
    {
        mSource = source;
        mSearchDistance = SvConstants.DISCORDANTPAIRSEARCHDISTANCE;
        mMinQuality = SvConstants.DISCORDANTPAIRMINMAPQ;
        mMinFragmentLength = SvConstants.DISCORDANTPAIRFRAGMENTLENGTH;
    }

    public List<Record> findDiscordantReads(final SizedIterable<Record> existingRecords)
    {
        if(!SvConstants.TRYEXTENDINGUSINGDISCORDANTREADS)
            return List.of();

        final Set<String> excludeFragments = existingRecords.stream()
                .map(Record::getName)
                .collect(Collectors.toSet());

        final Stream<RegionOfInterest> regionsUnmerged = existingRecords.stream().flatMap(this::interestingRegions);
        final List<RegionOfInterest> mappedSearchRegions = RegionOfInterest.tryMerge(regionsUnmerged::iterator);

        final List<Record> discordantReads = mappedSearchRegions.stream()
                .flatMap(region -> mSource.findReadsContaining(region.Chromosome, region.Start, region.End).stream())
                .filter(record -> !excludeFragments.contains(record.getName()))
                .filter(record -> record.isDiscordant(mMinFragmentLength))
                .collect(Collectors.toList());

        final List<Record> filteredDiscordantReads = filterFewMateRegions(filterHighDepthRegions(discordantReads));
        return filteredDiscordantReads.stream().distinct().collect(Collectors.toList());
    }

    private List<Record> filterFewMateRegions(final List<Record> records)
    {
        if(records.size() < 20)
            return records;

        final List<Pair<RegionOfInterest, List<Record>>> regionsWithDepth = tryMergeMateRegions(records);
        while(regionsWithDepth.size() > 8)
        {
            final int worstDepth = regionsWithDepth.stream()
                    .mapToInt(pair -> pair.getRight().size())
                    .min().orElseThrow();
            regionsWithDepth.removeIf(pair -> pair.getRight().size() <= worstDepth);
        }

        return regionsWithDepth.stream()
                .flatMap(pair -> pair.getRight().stream())
                .collect(Collectors.toList());
    }

    private List<Record> filterHighDepthRegions(final List<Record> records)
    {
        if(records.size() < 1_000)
            return records;

        final List<Pair<RegionOfInterest, List<Record>>> regionsWithDepth = tryMergeMateRegions(records);
        return regionsWithDepth.stream()
                .filter(pair -> pair.getRight().size() < 1_000)
                .flatMap(pair -> pair.getRight().stream())
                .collect(Collectors.toList());
    }

    private static List<Pair<RegionOfInterest, List<Record>>> tryMergeMateRegions(final List<Record> records)
    {
        final Map<String, List<Record>> byChromosome = new HashMap<>();
        for(final Record record : records)
        {
            if(!record.isMateMapped())
                continue;

            byChromosome.computeIfAbsent(record.getMateChromosome(), ignored -> new ArrayList<>())
                    .add(record);
        }

        // Sorting ensures that only a single pass is required to combine the regions, and that we only have to check the last one
        byChromosome.values().forEach(list -> list.sort(Comparator.comparingInt(Record::getMateAlignmentStart)));

        final Map<RegionOfInterest, List<Record>> depth = new IdentityHashMap<>();
        for(final List<Record> recordList : byChromosome.values())
        {
            final List<RegionOfInterest> consolidated = new ArrayList<>();
            for(final Record record : recordList)
            {
                final RegionOfInterest region = new RegionOfInterest(record.getMateChromosome(),
                        record.getMateAlignmentStart(), record.getMateAlignmentStart() + record.getLength());

                boolean found = false;
                for(final RegionOfInterest existingRegion : consolidated)
                    if(existingRegion.tryExtend(region))
                    {
                        depth.computeIfAbsent(existingRegion, __ -> new ArrayList<>()).add(record);
                        found = true;
                        break;
                    }

                if(!found)
                {
                    depth.computeIfAbsent(region, __ -> new ArrayList<>()).add(record);
                    consolidated.add(region);
                }
            }
        }

        return depth.entrySet().stream()
                .map(entry -> Pair.of(entry.getKey(), entry.getValue()))
                .collect(Collectors.toList());
    }

    private Stream<RegionOfInterest> interestingRegions(final Record record)
    {
        if(record.getMappingQuality() < mMinQuality)
            return Stream.of();

        final RegionOfInterest primary = new RegionOfInterest(record.getChromosome(),
                record.getUnclippedStart() - mSearchDistance, record.getUnclippedEnd() + mSearchDistance);
        if(record.isMateMapped())
        {
            final RegionOfInterest secondary = new RegionOfInterest(
                    record.getMateChromosome(),
                    record.getMateAlignmentStart() - mSearchDistance,
                    record.getMateAlignmentStart() + record.getLength() + mSearchDistance);
            return Stream.of(primary, secondary);
        }
        return Stream.of(primary);
    }
}
