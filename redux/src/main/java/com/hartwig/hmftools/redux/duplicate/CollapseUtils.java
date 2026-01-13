package com.hartwig.hmftools.redux.duplicate;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.redux.duplicate.FragmentCoords.FRAG_ORIENT_FORWARD_STR;
import static com.hartwig.hmftools.redux.duplicate.FragmentCoords.FRAG_ORIENT_REVERSE_STR;
import static com.hartwig.hmftools.redux.duplicate.FragmentCoords.FRAG_READ_LOWER_STR;
import static com.hartwig.hmftools.redux.duplicate.FragmentCoords.FRAG_READ_UPPER_STR;
import static com.hartwig.hmftools.redux.duplicate.FragmentCoords.FRAG_TYPE_SUPP_INFO_STR;
import static com.hartwig.hmftools.redux.duplicate.FragmentCoords.FRAG_TYPE_UNMAPPED;

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
                singleReads.add(new ReadInfo(collapsedGroup.reads().get(0), collapsedGroup.fragmentCoordinates()));
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

        String lowerOrientation = fragmentCoords.OrientLower == FORWARD ? FRAG_ORIENT_FORWARD_STR : FRAG_ORIENT_REVERSE_STR;
        String suppSuffix = fragmentCoords.SuppReadInfo == null ? "" : ":" + FRAG_TYPE_SUPP_INFO_STR;
        if(fragmentCoords.PositionUpper == NO_POSITION)
        {
            String unmappedSuffix = fragmentCoords.UnmappedSourced ? format(":%c", FRAG_TYPE_UNMAPPED) : "";
            return format("%s:%s%s%s", fragmentCoords.ChromsomeLower, lowerOrientation, unmappedSuffix, suppSuffix);
        }

        String upperOrientation = fragmentCoords.OrientUpper == FORWARD ? FRAG_ORIENT_FORWARD_STR : FRAG_ORIENT_REVERSE_STR;
        String isLowerString = fragmentCoords.ReadIsLower ? FRAG_READ_LOWER_STR : FRAG_READ_UPPER_STR;
        return format("%s:%s:%s:%s:%s%s",
                fragmentCoords.ChromsomeLower, lowerOrientation, fragmentCoords.ChromsomeUpper, upperOrientation, isLowerString, suppSuffix);
    }
}
