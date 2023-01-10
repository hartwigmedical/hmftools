package com.hartwig.hmftools.bamtools.markdups;

import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

public class UmiGroupUtils
{
    public static final int MAX_UMI_BASE_DIFF = 1;

    private static final String READ_ID_DELIM = ":";

    public static boolean exceedsUmiIdDiff(final String first, final String second)
    {
        if(first.length() != second.length())
            return true;

        short diffs = 0;
        for(short i = 0; i < first.length(); ++i)
        {
            if(first.charAt(i) != second.charAt(i))
            {
                ++diffs;

                if(diffs > MAX_UMI_BASE_DIFF)
                    return true;
            }
        }

        return false;
    }

    public static String extractUmiId(final String readId)
    {
        String[] items = readId.split(READ_ID_DELIM, -1);
        return items[items.length - 1];
    }

    private class UmiGroupCluster
    {
        public final List<UmiGroup> Groups;

        public UmiGroupCluster(final UmiGroup group)
        {
            Groups = Lists.newArrayList(group);
        }
    }

    public static List<UmiGroup> buildUmiGroups(final List<Fragment> fragments)
    {
        Map<String,UmiGroup> groups = Maps.newHashMap();

        for(Fragment fragment : fragments)
        {
            String umiId = extractUmiId(fragment.id());

            UmiGroup group = groups.get(umiId);

            if(group == null)
            {
                groups.put(umiId, new UmiGroup(umiId, fragment));
            }
            else
            {
                group.Fragments.add(fragment);
            }
        }

        List<UmiGroup> orderedGroups = groups.values().stream().sorted(new SizeComparator()).collect(Collectors.toList());
        List<UmiGroup> finalGroups = Lists.newArrayList();

        int i = 0;
        while(i < orderedGroups.size() - 1)
        {
            UmiGroup first = orderedGroups.get(i);

            List<UmiGroup> cluster = Lists.newArrayList(first);

            int j = i + 1;
            while(j < orderedGroups.size())
            {
                UmiGroup second = orderedGroups.get(j);

                boolean merged = false;

                for(UmiGroup existing : cluster)
                {
                    if(existing.fragmentCount() > second.fragmentCount() && !exceedsUmiIdDiff(existing.UmiId, second.UmiId))
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
                }
            }

            for(j = 1; j < cluster.size(); ++j)
            {
                first.Fragments.addAll(cluster.get(j).Fragments);
            }

            finalGroups.add(first);
            ++i;
        }

        return orderedGroups;
    }

    private static class SizeComparator implements Comparator<UmiGroup>
    {
        public int compare(final UmiGroup first, final UmiGroup second)
        {
            return first.Fragments.size() < second.Fragments.size() ? 1 : -1;
        }
    }


}
