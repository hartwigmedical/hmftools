package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.esvee.read.ReadUtils.isDiscordant;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.esvee.common.RegionOfInterest;
import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.read.SAMSource;
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

    public List<Read> findDiscordantReads(final SizedIterable<Read> existingRecords)
    {
        if(!SvConstants.TRYEXTENDINGUSINGDISCORDANTREADS)
            return List.of();

        final Set<String> excludeFragments = existingRecords.stream()
                .map(Read::getName)
                .collect(Collectors.toSet());

        final Stream<RegionOfInterest> regionsUnmerged = existingRecords.stream().flatMap(this::interestingRegions);
        final List<RegionOfInterest> mappedSearchRegions = RegionOfInterest.tryMerge(regionsUnmerged::iterator);

        final List<Read> discordantReads = mappedSearchRegions.stream()
                .flatMap(region -> mSource.findReadsContaining(region.Chromosome, region.Start, region.End).stream())
                .filter(record -> !excludeFragments.contains(record.getName()))
                .filter(record -> isDiscordant(record, mMinFragmentLength))
                .collect(Collectors.toList());

        final List<Read> filteredDiscordantReads = filterFewMateRegions(filterHighDepthRegions(discordantReads));
        return filteredDiscordantReads.stream().distinct().collect(Collectors.toList());
    }

    private List<Read> filterFewMateRegions(final List<Read> reads)
    {
        if(reads.size() < 20)
            return reads;

        final List<Pair<RegionOfInterest, List<Read>>> regionsWithDepth = tryMergeMateRegions(reads);
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

    private List<Read> filterHighDepthRegions(final List<Read> reads)
    {
        if(reads.size() < 1_000)
            return reads;

        final List<Pair<RegionOfInterest, List<Read>>> regionsWithDepth = tryMergeMateRegions(reads);
        return regionsWithDepth.stream()
                .filter(pair -> pair.getRight().size() < 1_000)
                .flatMap(pair -> pair.getRight().stream())
                .collect(Collectors.toList());
    }

    private static List<Pair<RegionOfInterest, List<Read>>> tryMergeMateRegions(final List<Read> reads)
    {
        final Map<String, List<Read>> byChromosome = new HashMap<>();
        for(Read read : reads)
        {
            if(!read.isMateMapped())
                continue;

            byChromosome.computeIfAbsent(read.getMateChromosome(), ignored -> new ArrayList<>())
                    .add(read);
        }

        // Sorting ensures that only a single pass is required to combine the regions, and that we only have to check the last one
        byChromosome.values().forEach(list -> list.sort(Comparator.comparingInt(Read::getMateAlignmentStart)));

        final Map<RegionOfInterest, List<Read>> depth = new IdentityHashMap<>();
        for(List<Read> readList : byChromosome.values())
        {
            final List<RegionOfInterest> consolidated = new ArrayList<>();
            for(Read read : readList)
            {
                final RegionOfInterest region = new RegionOfInterest(read.getMateChromosome(),
                        read.getMateAlignmentStart(), read.getMateAlignmentStart() + read.getLength());

                boolean found = false;
                for(RegionOfInterest existingRegion : consolidated)
                    if(existingRegion.tryExtend(region))
                    {
                        depth.computeIfAbsent(existingRegion, __ -> new ArrayList<>()).add(read);
                        found = true;
                        break;
                    }

                if(!found)
                {
                    depth.computeIfAbsent(region, __ -> new ArrayList<>()).add(read);
                    consolidated.add(region);
                }
            }
        }

        return depth.entrySet().stream()
                .map(entry -> Pair.of(entry.getKey(), entry.getValue()))
                .collect(Collectors.toList());
    }

    private Stream<RegionOfInterest> interestingRegions(final Read read)
    {
        if(read.getMappingQuality() < mMinQuality)
            return Stream.of();

        final RegionOfInterest primary = new RegionOfInterest(read.getChromosome(),
                read.getUnclippedStart() - mSearchDistance, read.getUnclippedEnd() + mSearchDistance);
        if(read.isMateMapped())
        {
            final RegionOfInterest secondary = new RegionOfInterest(
                    read.getMateChromosome(),
                    read.getMateAlignmentStart() - mSearchDistance,
                    read.getMateAlignmentStart() + read.getLength() + mSearchDistance);
            return Stream.of(primary, secondary);
        }
        return Stream.of(primary);
    }
}
