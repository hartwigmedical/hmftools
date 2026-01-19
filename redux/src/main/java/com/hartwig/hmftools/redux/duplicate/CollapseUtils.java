package com.hartwig.hmftools.redux.duplicate;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.redux.duplicate.FragmentCoords.COORD_ORIENT_FORWARD_STR;
import static com.hartwig.hmftools.redux.duplicate.FragmentCoords.COORD_ORIENT_REVERSE_STR;
import static com.hartwig.hmftools.redux.duplicate.FragmentCoords.COORD_READ_LOWER_STR;
import static com.hartwig.hmftools.redux.duplicate.FragmentCoords.COORD_READ_UPPER_STR;
import static com.hartwig.hmftools.redux.duplicate.FragmentCoords.COORD_READ_SUPP_INFO_STR;
import static com.hartwig.hmftools.redux.duplicate.FragmentCoords.COORD_READ_UNMAPPED;

import java.util.Collection;
import java.util.List;
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
                singleReads.add(new ReadInfo(collapsedGroup.reads().get(0), collapsedGroup.fragCoordinates()));
                return;
            }

            duplicateGroups.add(collapsedGroup);
        });

        return new FragmentCoordReads(duplicateGroups, singleReads);
    }

    public static FragmentCoordReads getFragmentCoordReads(final Collection<DuplicateGroup> collapsedGroups)
    {
        return getFragmentCoordReads(collapsedGroups.stream());
    }

    @Nullable
    static String collapseToNonOrientedKeyWithoutCoordinates(final FragmentCoords fragmentCoords)
    {
        if(fragmentCoords.Unpaired)
            return null;

        String lowerOrientation = fragmentCoords.OrientLower == FORWARD ? COORD_ORIENT_FORWARD_STR : COORD_ORIENT_REVERSE_STR;
        String suppSuffix = fragmentCoords.SuppReadInfo == null ? "" : ":" + COORD_READ_SUPP_INFO_STR;
        if(fragmentCoords.PositionUpper == NO_POSITION)
        {
            String unmappedSuffix = fragmentCoords.UnmappedSourced ? format(":%c", COORD_READ_UNMAPPED) : "";
            return format("%s:%s%s%s", fragmentCoords.ChromsomeLower, lowerOrientation, unmappedSuffix, suppSuffix);
        }

        String upperOrientation = fragmentCoords.OrientUpper == FORWARD ? COORD_ORIENT_FORWARD_STR : COORD_ORIENT_REVERSE_STR;
        String isLowerString = fragmentCoords.ReadIsLower ? COORD_READ_LOWER_STR : COORD_READ_UPPER_STR;
        return format("%s:%s:%s:%s:%s%s",
                fragmentCoords.ChromsomeLower, lowerOrientation, fragmentCoords.ChromsomeUpper, upperOrientation, isLowerString, suppSuffix);
    }
}
