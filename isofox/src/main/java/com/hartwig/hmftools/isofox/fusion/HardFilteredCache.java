package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.chromosomeRank;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.fusion.SupplementaryJunctionData.fromReads;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class HardFilteredCache
{
    private final Map<String,Set<String>> mChromosomePairFilteredReads;

    public HardFilteredCache()
    {
        mChromosomePairFilteredReads = Maps.newHashMap();
    }

    public int cacheCount() { return mChromosomePairFilteredReads.values().stream().mapToInt(x -> x.size()).sum(); }
    public int chrPairCount() { return mChromosomePairFilteredReads.size(); }

    public void addHardFilteredReads(final Map<String,Set<String>> hardFilteredReadIds)
    {
        for(Map.Entry<String,Set<String>> entry : hardFilteredReadIds.entrySet())
        {
            if(entry.getValue().isEmpty())
                continue;

            Set<String> filteredReadIds = mChromosomePairFilteredReads.get(entry.getKey());
            if(filteredReadIds == null)
            {
                filteredReadIds = Sets.newHashSet();
                mChromosomePairFilteredReads.put(entry.getKey(), filteredReadIds);
            }

            filteredReadIds.addAll(entry.getValue());
        }
    }

    public void removeHardFilteredReads(
            final String chromosome, final Map<String,Map<String,FusionReadGroup>> chrIncompleteGroups,
            final Map<String,Set<String>> chrHardFilteredReadIds)
    {
        // removes partial groups which were hard-filtered by other chromosomes
        // and removes these same reads from the new hard-filtered groups
        // and reconciles matching hard-filtered reads from each chromosome (the expected the scenario)
        for(Map.Entry<String,Map<String,FusionReadGroup>> chrEntry : chrIncompleteGroups.entrySet())
        {
            String otherChromosome = chrEntry.getKey();
            String chrPair = formChromosomePairString(chromosome, otherChromosome);

            Map<String,FusionReadGroup> incompleteGroups = chrEntry.getValue();

            List<String> incompleteFilteredReadIds = incompleteGroups.keySet().stream()
                    .filter(x -> wasHardFiltered(chrPair, x)).collect(Collectors.toList());

            incompleteFilteredReadIds.forEach(x -> incompleteGroups.remove(x));

            Set<String> newHardFilteredReads = chrHardFilteredReadIds.get(chrPair);
            Set<String> existingFilteredReadIds = mChromosomePairFilteredReads.get(chrPair);

            // reconcile filtered reads, expecting these to find a match
            if(newHardFilteredReads != null && existingFilteredReadIds != null)
            {
                Set<String> matchedReadIds = newHardFilteredReads.stream().filter(x -> existingFilteredReadIds.contains(x)).collect(Collectors.toSet());
                matchedReadIds.forEach(x -> newHardFilteredReads.remove(x));
                matchedReadIds.forEach(x -> existingFilteredReadIds.remove(x));
            }

            // remove new filtered groups if they match a previously hard-filtered read id
            if(newHardFilteredReads != null)
                incompleteFilteredReadIds.forEach(x -> newHardFilteredReads.remove(x));
        }
    }

    public void purgeChromosomeEntries(final String chromosome, final String otherChromosome)
    {
        String chrPair = formChromosomePairString(chromosome, otherChromosome);
        mChromosomePairFilteredReads.remove(chrPair);
    }

    public static void removePartialGroupsWithHardFilteredMatch(
            final Map<String,FusionReadGroup> partialGroups, final Map<String,Set<String>> newHardFilteredReadIds)
    {
        // remove any partial groups with a hard-filter match, and vice versa
        for(Set<String> hardFilteredReads : newHardFilteredReadIds.values())
        {
            Set<String> matchedReadIds = hardFilteredReads.stream().filter(x -> partialGroups.containsKey(x)).collect(Collectors.toSet());
            matchedReadIds.forEach(x -> hardFilteredReads.remove(x));
            matchedReadIds.forEach(x -> partialGroups.remove(x));
        }
    }

    private static boolean containsChromosome(final String chrPair, final String chromosome)
    {
        String[] chromosomes = chrPair.split("_", 2);
        return chromosomes[0].equals(chromosome) || chromosomes[1].equals(chromosome);
    }

    private boolean wasHardFiltered(final String chrPair, final String readId)
    {
        Set<String> filteredReadIds = mChromosomePairFilteredReads.get(chrPair);
        return filteredReadIds != null && filteredReadIds.contains(readId);
    }

    public static String formChromosomePairString(final String chr1, final String chr2)
    {
        int chrRank1 = chromosomeRank(chr1);
        int chrRank2 = chromosomeRank(chr2);
        String[] chromosomes = new String[SE_PAIR];

        if(chrRank1 <= chrRank2)
        {
            chromosomes[SE_START] = chr1;
            chromosomes[SE_END] = chr2;
        }
        else
        {
            chromosomes[SE_START] = chr2;
            chromosomes[SE_END] = chr1;
        }

        return FusionUtils.formChromosomePair(chromosomes);
    }

    public static void applyHardFilter(
            final Map<Integer,List<SupplementaryJunctionData>> supplementaryJunctions, final Map<String,Set<String>> hardFilteredReadIds,
            final Map<String,ChimericReadGroup> chimericReadMap, final String localChromosome, int minSplitFrags)
    {
        if(minSplitFrags <= 1)
            return;

        // filter will be repeated in the fusion finding routine when all info is known, but for now hard-filter any read groups
        // without a read matching a known fusion gene

        for(List<SupplementaryJunctionData> suppJunctions : supplementaryJunctions.values())
        {
            for(SupplementaryJunctionData suppJuncData : suppJunctions)
            {
                if(suppJuncData.MatchCount + 1 >= minSplitFrags)
                    continue;

                String chrPair = formChromosomePairString(localChromosome, suppJuncData.RemoteChromosome);

                Set<String> filteredReadIds = hardFilteredReadIds.get(chrPair);
                if(filteredReadIds == null)
                {
                    filteredReadIds = com.beust.jcommander.internal.Sets.newHashSet();
                    hardFilteredReadIds.put(chrPair, filteredReadIds);
                }

                for(String readId : suppJuncData.ReadIds)
                {
                    filteredReadIds.add(readId);
                    chimericReadMap.remove(readId);
                }
            }
        }
    }
}
