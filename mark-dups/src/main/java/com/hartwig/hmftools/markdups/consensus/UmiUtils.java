package com.hartwig.hmftools.markdups.consensus;

import static com.hartwig.hmftools.markdups.common.Constants.DEFAULT_MAX_UMI_BASE_DIFF;
import static com.hartwig.hmftools.markdups.common.FragmentCoordinates.FRAGMENT_REVERSED_ID;

import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.markdups.common.Fragment;

public final class UmiUtils
{
    public static List<DuplicateGroup> buildUmiGroups(final List<Fragment> fragments, final UmiConfig config)
    {
        Map<String, DuplicateGroup> groups = Maps.newHashMap();
        boolean checkDefinedUmis = config.hasDefinedUmis();
        boolean useDefinedUmis = checkDefinedUmis;

        for(Fragment fragment : fragments)
        {
            String umiId = config.extractUmiId(fragment.id());

            if(checkDefinedUmis)
            {
                String definedUmiId = config.matchDefinedUmiId(umiId);
                if(definedUmiId == null)
                {
                    useDefinedUmis = false;
                    checkDefinedUmis = false;
                }
                else
                {
                    umiId = definedUmiId;
                }
            }

            DuplicateGroup group = groups.get(umiId);

            if(group == null)
            {
                groups.put(umiId, new DuplicateGroup(umiId, fragment));
            }
            else
            {
                group.fragments().add(fragment);
            }
        }

        if(useDefinedUmis)
        {
            return groups.values().stream().collect(Collectors.toList());
        }

        // order groups by descending number of fragments
        List<DuplicateGroup> orderedGroups = groups.values().stream().sorted(new SizeComparator()).collect(Collectors.toList());

        // then apply the directional model, where smaller groups are merged into larger ones
        int i = 0;
        while(i < orderedGroups.size() - 1)
        {
            DuplicateGroup first = orderedGroups.get(i);

            List<DuplicateGroup> cluster = Lists.newArrayList(first);

            int j = i + 1;
            while(j < orderedGroups.size())
            {
                DuplicateGroup second = orderedGroups.get(j);

                boolean merged = false;

                for(DuplicateGroup existing : cluster)
                {
                    if(existing.fragmentCount() >= second.fragmentCount() && !exceedsUmiIdDiff(existing.id(), second.id(), config.PermittedBaseDiff))
                    {
                        merged = true;
                        break;
                    }
                }

                if(!merged)
                {
                    ++j;
                }
                else
                {
                    orderedGroups.remove(j);
                    cluster.add(second);

                    // restart the search since a newly added group may be close enough to a skipped one
                    j = i + 1;
                }
            }

            for(j = 1; j < cluster.size(); ++j)
            {
                first.fragments().addAll(cluster.get(j).fragments());
            }

            ++i;
        }

        // run a final check allowing collapsing of UMIs with 2-base differences
        if(orderedGroups.size() > 1)
        {
            i = 0;
            while(i < orderedGroups.size())
            {
                DuplicateGroup first = orderedGroups.get(i);

                int j = i + 1;
                while(j < orderedGroups.size())
                {
                    DuplicateGroup second = orderedGroups.get(j);

                    if(!exceedsUmiIdDiff(first.id(), second.id(), config.PermittedBaseDiff + 1))
                    {
                        first.fragments().addAll(second.fragments());
                        orderedGroups.remove(j);
                    }
                    else
                    {
                        ++j;
                    }
                }

                ++i;
            }
        }

        return orderedGroups;
    }

    public static boolean exceedsUmiIdDiff(final String first, final String second)
    {
        return exceedsUmiIdDiff(first, second, DEFAULT_MAX_UMI_BASE_DIFF);
    }

    public static boolean exceedsUmiIdDiff(final String first, final String second, int permittedDiff)
    {
        if(first.length() != second.length())
            return true;

        short diffs = 0;
        for(short i = 0; i < first.length(); ++i)
        {
            if(first.charAt(i) != second.charAt(i))
            {
                ++diffs;

                if(diffs > permittedDiff)
                    return true;
            }
        }

        return false;
    }

    public static int calcUmiIdDiff(final String first, final String second)
    {
        if(first.length() != second.length())
            return 0;

        int diffs = 0;
        for(short i = 0; i < first.length(); ++i)
        {
            if(first.charAt(i) != second.charAt(i))
            {
                ++diffs;
            }
        }

        return diffs;
    }

    public static void collapseOnDuplexMatches(final List<DuplicateGroup> duplicateGroups, final UmiConfig config)
    {
        int i = 0;
        while(i < duplicateGroups.size())
        {
            DuplicateGroup first = duplicateGroups.get(i);

            int j = i + 1;
            while(j < duplicateGroups.size())
            {
                DuplicateGroup second = duplicateGroups.get(j);

                if(hasReversedFragmentCoords(first.coordinatesKey(), second.coordinatesKey())
                && hasDuplexUmiMatch(first.id(), second.id(), config.DuplexDelim, config.PermittedBaseDiff))
                {
                    first.fragments().addAll(second.fragments());
                    first.registerDuplex();
                    duplicateGroups.remove(j);
                }
                else
                {
                    ++j;
                }
            }

            ++i;
        }
    }

    public static boolean hasReversedFragmentCoords(final String coords1, final String coords2)
    {
        int coord1ReversedIndex = coords1.indexOf(FRAGMENT_REVERSED_ID);
        int coord2ReversedIndex = coords2.indexOf(FRAGMENT_REVERSED_ID);

        if(coord1ReversedIndex == coord2ReversedIndex)
            return false;

        if(coord1ReversedIndex > 0)
            return coords1.substring(0, coord1ReversedIndex - 1).equals(coords2);
        else
            return coords2.substring(0, coord2ReversedIndex - 1).equals(coords1);
    }

    public static boolean hasDuplexUmiMatch(final String first, final String second, final String duplexDelim, int permittedDiff)
    {
        String[] umiParts1 = first.split(duplexDelim, 2);
        String[] umiParts2 = second.split(duplexDelim, 2);

        if(umiParts1.length != 2 || umiParts2.length != 2)
            return false;

        return !exceedsUmiIdDiff(umiParts1[0], umiParts2[1], permittedDiff) && !exceedsUmiIdDiff(umiParts1[1], umiParts2[0], permittedDiff);
    }

    private static class SizeComparator implements Comparator<DuplicateGroup>
    {
        public int compare(final DuplicateGroup first, final DuplicateGroup second)
        {
            if(first.fragments().size() < second.fragments().size())
                return 1;
            else if(first.fragments().size() > second.fragments().size())
                return -1;
            else
                return 0;
        }
    }

}
