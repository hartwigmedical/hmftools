package com.hartwig.hmftools.bamtools.markdups;

import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.PRIMARY;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.findPrimaryFragment;

import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class DuplicateGroupUtils
{
    public static final int MAX_UMI_BASE_DIFF = 1;

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

    public static void processDuplicateGroups(final List<List<Fragment>> duplicateGroups, final UmiConfig umiConfig)
    {
        if(duplicateGroups == null)
            return;

        if(umiConfig.Enabled)
        {
            for(List<Fragment> fragments : duplicateGroups)
            {
                List<UmiGroup> umiGroups = buildUmiGroups(fragments, umiConfig);

                for(UmiGroup umiGroup : umiGroups)
                {
                    umiGroup.Fragments.forEach(x -> x.setUmiId(umiGroup.UmiId));

                    setDuplicateGroupProperties(umiGroup.Fragments);
                }
            }
        }
        else
        {
            for(List<Fragment> fragments : duplicateGroups)
            {
                setDuplicateGroupProperties(fragments);
            }
        }
    }

    private static void setDuplicateGroupProperties(final List<Fragment> duplicateFragments)
    {
        int dupCount = duplicateFragments.size();
        duplicateFragments.forEach(x -> x.setDuplicateCount(dupCount));

        Fragment primary = findPrimaryFragment(duplicateFragments, false);
        primary.setStatus(PRIMARY);
    }


    public static List<UmiGroup> buildUmiGroups(final List<Fragment> fragments, final UmiConfig config)
    {
        Map<String,UmiGroup> groups = Maps.newHashMap();

        for(Fragment fragment : fragments)
        {
            String umiId = config.extractUmiId(fragment.id());

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
