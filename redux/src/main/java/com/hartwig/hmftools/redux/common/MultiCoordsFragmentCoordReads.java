package com.hartwig.hmftools.redux.common;

import static java.lang.String.format;

import java.util.Collection;
import java.util.List;

import com.google.common.collect.Lists;

public class MultiCoordsFragmentCoordReads extends FragmentCoordReads
{
    public final List<MultiCoordsDuplicateGroup> MultiCoordDuplicateGroups;

    public MultiCoordsFragmentCoordReads(final List<DuplicateGroup> duplicateGroups, final List<ReadInfo> singleReads,
            final List<MultiCoordsDuplicateGroup> multiCoordDuplicateGroups)
    {
        super(duplicateGroups, singleReads);

        MultiCoordDuplicateGroups = multiCoordDuplicateGroups != null ? multiCoordDuplicateGroups : Lists.newArrayList();
    }

    public int multiCoordDuplicateGroupReadCount()
    {
        return MultiCoordDuplicateGroups.stream().mapToInt(x -> x.readCount()).sum();
    }

    @Override
    public int totalReadCount()
    {
        return super.totalReadCount() + multiCoordDuplicateGroupReadCount();
    }

    @Override
    public int coordinateCount()
    {
        return super.coordinateCount() + MultiCoordDuplicateGroups.size();
    }

    @Override
    public int minReadPositionStart()
    {
        int minPositionStart = super.minReadPositionStart();

        for(MultiCoordsDuplicateGroup duplicateGroup : MultiCoordDuplicateGroups)
        {
            int minGroupReadPosition = duplicateGroup.reads().stream().mapToInt(x -> x.getAlignmentStart()).min().orElse(-1);

            if(minPositionStart < 0 || minGroupReadPosition > 0 && minGroupReadPosition < minPositionStart)
                minPositionStart = minGroupReadPosition;
        }

        return minPositionStart;
    }

    public static FragmentCoordReads fromCollapsedGroups(final Collection<MultiCoordsDuplicateGroup> collapsedGroups)
    {
        List<ReadInfo> singleReads = Lists.newArrayList();
        List<DuplicateGroup> duplicateGroups = Lists.newArrayList();
        List<MultiCoordsDuplicateGroup> multiCoordGroups = Lists.newArrayList();

        for(MultiCoordsDuplicateGroup collapsedGroup : collapsedGroups)
        {
            if(collapsedGroup.isSingleRead())
            {
                singleReads.add(collapsedGroup.toReadInfo());
                continue;
            }

            if(collapsedGroup.isSingleCoordDuplicateGroup())
            {
                duplicateGroups.add(collapsedGroup.toSingleCoordDuplicateGroup());
                continue;
            }

            multiCoordGroups.add(collapsedGroup);
        }

        if(multiCoordGroups.isEmpty())
            return new FragmentCoordReads(duplicateGroups, singleReads);

        return new MultiCoordsFragmentCoordReads(duplicateGroups, singleReads, multiCoordGroups);
    }

    @Override
    public String toString()
    {
        return format("duplicateGroups(%d) singles(%d) multiCoordDuplicateGroups(%d) totalReads(%d)",
                DuplicateGroups.size(), SingleReads.size(), MultiCoordDuplicateGroups.size(), totalReadCount());
    }
}
