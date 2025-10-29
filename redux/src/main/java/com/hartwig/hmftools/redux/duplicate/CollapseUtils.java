package com.hartwig.hmftools.redux.duplicate;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;

import java.util.Collection;
import java.util.List;
import java.util.function.BinaryOperator;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.redux.common.ReadInfo;

import org.jetbrains.annotations.Nullable;

public final class CollapseUtils
{
    public static FragmentCoordReads getFragmentCoordReads(final Stream<DuplicateGroup> collapsedGroups)
    {
        final List<DuplicateGroup> duplicateGroups = Lists.newArrayList();
        final List<ReadInfo> singleReads = Lists.newArrayList();
        collapsedGroups.forEach(collapsedGroup ->
        {
            if(collapsedGroup.totalReadCount() == 1)
            {
                singleReads.add(new ReadInfo(collapsedGroup.reads().get(0), collapsedGroup.fragmentCoordinates()));
                return;
            }

            duplicateGroups.add(collapsedGroup);
        });

        return new FragmentCoordReads(duplicateGroups, singleReads);
    }

    public static BinaryOperator<DuplicateGroup> DUPLICATE_GROUP_MERGER = (acc, group) ->
    {
        acc.addReads(group.reads());
        return acc;
    };

    public static FragmentCoordReads getFragmentCoordReads(final Collection<DuplicateGroup> collapsedGroups)
    {
        return getFragmentCoordReads(collapsedGroups.stream());
    }

    @Nullable
    static String collapseToNonOrientedKeyWithoutCoordinates(final FragmentCoords fragmentCoords)
    {
        if(fragmentCoords.Unpaired)
            return null;

        String lowerOrientation = fragmentCoords.OrientLower == FORWARD ? "F" : "R";
        String suppSuffix = fragmentCoords.SuppReadInfo == null ? "" : ":S";
        if(fragmentCoords.PositionUpper == NO_POSITION)
        {
            String unmappedSuffix = fragmentCoords.UnmappedSourced ? ":U" : "";
            return format("%s:%s%s%s", fragmentCoords.ChromsomeLower, lowerOrientation, unmappedSuffix, suppSuffix);
        }

        String upperOrientation = fragmentCoords.OrientUpper == FORWARD ? "F" : "R";
        String isLowerString = fragmentCoords.ReadIsLower ? "L" : "U";
        return format("%s:%s:%s:%s:%s%s",
                fragmentCoords.ChromsomeLower, lowerOrientation, fragmentCoords.ChromsomeUpper, upperOrientation, isLowerString, suppSuffix);
    }
}
