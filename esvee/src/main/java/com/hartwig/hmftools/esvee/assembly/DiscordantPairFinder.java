package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.esvee.read.ReadUtils.isDiscordant;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.esvee.old.RegionOfInterest;
import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.read.Read;

import org.apache.commons.lang3.tuple.Pair;

public class DiscordantPairFinder
{
    private final int mSearchDistance;
    private final int mMinQuality;
    private final int mMinFragmentLength;

    public DiscordantPairFinder()
    {
        mSearchDistance = SvConstants.DISCORDANT_PAIR_SEARCH_DISTANCE;
        mMinQuality = SvConstants.DISCORDANT_PAIR_MIN_MAPQ;
        mMinFragmentLength = SvConstants.DISCORDANT_FRAGMENT_LENGTH;
    }

    public List<Read> findDiscordantReads(final List<Read> existingSupportingReads, final List<Read> allReads)
    {
        if(!SvConstants.TRY_EXTENDING_USING_DISCORDANT_READS)
            return Collections.emptyList();

        List<Read> discordantReads = allReads.stream()
                .filter(x -> !existingSupportingReads.contains(x))
                .filter(x -> isDiscordant(x))
                .collect(Collectors.toList());

        List<Read> filteredDiscordantReads = filterFewMateRegions(filterHighDepthRegions(discordantReads));
        return filteredDiscordantReads.stream().distinct().collect(Collectors.toList());
    }

    /*
    public List<Read> findDiscordantReads(final SizedIterable<Read> existingRecords)
    {
        if(!SvConstants.TRY_EXTENDING_USING_DISCORDANT_READS)
            return List.of();

        final Set<String> excludeFragments = existingRecords.stream()
                .map(Read::getName)
                .collect(Collectors.toSet());

        final Stream<RegionOfInterest> regionsUnmerged = existingRecords.stream().flatMap(this::interestingRegions);
        final List<RegionOfInterest> mappedSearchRegions = RegionOfInterest.tryMerge(regionsUnmerged::iterator);

        final List<Read> discordantReads = mappedSearchRegions.stream()
                .flatMap(region -> mSource.findReadsContaining(region.Chromosome, region.start(), region.end()).stream())
                .filter(record -> !excludeFragments.contains(record.getName()))
                .filter(record -> isDiscordant(record, mMinFragmentLength))
                .collect(Collectors.toList());

        final List<Read> filteredDiscordantReads = filterFewMateRegions(filterHighDepthRegions(discordantReads));
        return filteredDiscordantReads.stream().distinct().collect(Collectors.toList());
    }
    */

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

            byChromosome.computeIfAbsent(read.mateChromosome(), ignored -> new ArrayList<>())
                    .add(read);
        }

        // Sorting ensures that only a single pass is required to combine the regions, and that we only have to check the last one
        byChromosome.values().forEach(list -> list.sort(Comparator.comparingInt(Read::mateAlignmentStart)));

        final Map<RegionOfInterest, List<Read>> depth = new IdentityHashMap<>();
        for(List<Read> readList : byChromosome.values())
        {
            final List<RegionOfInterest> consolidated = new ArrayList<>();
            for(Read read : readList)
            {
                final RegionOfInterest region = new RegionOfInterest(read.mateChromosome(),
                        read.mateAlignmentStart(), read.mateAlignmentStart() + read.basesLength());

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
        if(read.mappingQuality() < mMinQuality)
            return Stream.of();

        final RegionOfInterest primary = new RegionOfInterest(read.chromosome(),
                read.unclippedStart() - mSearchDistance, read.unclippedEnd() + mSearchDistance);
        if(read.isMateMapped())
        {
            final RegionOfInterest secondary = new RegionOfInterest(
                    read.mateChromosome(),
                    read.mateAlignmentStart() - mSearchDistance,
                    read.mateAlignmentStart() + read.basesLength() + mSearchDistance);
            return Stream.of(primary, secondary);
        }
        return Stream.of(primary);
    }
}
