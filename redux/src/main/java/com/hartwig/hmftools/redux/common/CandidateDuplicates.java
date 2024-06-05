package com.hartwig.hmftools.redux.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.redux.common.FragmentCoordinates.formCoordinate;
import static com.hartwig.hmftools.redux.common.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.redux.common.FragmentStatus.NONE;
import static com.hartwig.hmftools.redux.common.FragmentUtils.calcFragmentStatus;

import java.util.List;

import com.google.common.collect.Lists;

import htsjdk.samtools.SAMRecord;

public class CandidateDuplicates
{
    // incomplete fragments (ie missing a mate read) with a matching fragment coordinate, and so candidates for being duplicates
    private final String mKey;
    private final List<Fragment> mFragments;

    private boolean mFinalised;

    public CandidateDuplicates(final String key, final Fragment fragment)
    {
        mKey = key;
        mFragments = Lists.newArrayList();
        addFragment(fragment);
        mFinalised = false;
    }

    public static CandidateDuplicates from(final Fragment fragment)
    {
        final SAMRecord read = fragment.reads().get(0);
        boolean mateForwardStrand = !read.getMateNegativeStrandFlag();

        String key = format("%d_%s",
                fragment.initialPosition(), formCoordinate(read.getMateReferenceName(), read.getMateAlignmentStart(), mateForwardStrand));

        return new CandidateDuplicates(key, fragment);
    }

    public String key() { return mKey; }

    public List<Fragment> fragments() { return mFragments; }
    public int fragmentCount() { return mFragments.size(); }

    public void addFragment(final Fragment fragment)
    {
        mFragments.add(fragment);
        fragment.setCandidateDupKey(mKey);
    }

    public boolean allFragmentsReady() { return mFragments.stream().allMatch(x -> x.primaryReadsPresent()); }
    public boolean finalised() { return mFinalised; }

    public List<List<Fragment>> finaliseFragmentStatus(boolean requireOrientationMatch)
    {
        if(mFinalised || !allFragmentsReady())
            return null;

        mFinalised = true;

        if(mFragments.size() == 1)
        {
            Fragment fragment = mFragments.get(0);
            fragment.setStatus(NONE);
            return null;
        }

        List<List<Fragment>> duplicateGroups = null;

        for(int i = 0; i < mFragments.size(); ++i)
        {
            Fragment fragment1 = mFragments.get(i);

            if(fragment1.status().isDuplicate()) // already a part of a group
                continue;

            if(i == mFragments.size() - 1)
            {
                fragment1.setStatus(NONE);
                break;
            }

            List<Fragment> duplicateFragments = null;

            for(int j = i + 1; j < mFragments.size(); ++j)
            {
                Fragment fragment2 = mFragments.get(j);

                if(fragment2.status().isDuplicate()) // already a part of a group
                    continue;

                FragmentStatus status = calcFragmentStatus(fragment1, fragment2, requireOrientationMatch);

                if(status == DUPLICATE)
                {
                    fragment1.setStatus(status);
                    fragment2.setStatus(status);

                    if(duplicateFragments == null)
                        duplicateFragments = Lists.newArrayList(fragment1);

                    duplicateFragments.add(fragment2);
                }
            }

            if(fragment1.status().isDuplicate())
            {
                if(duplicateGroups == null)
                    duplicateGroups = Lists.newArrayList();

                duplicateGroups.add(duplicateFragments);
            }
            else
            {
                fragment1.setStatus(NONE);
            }
        }

        return duplicateGroups;
    }

    public String toString()
    {
        return format("%s: fragments(%d) finalised(%s)", mKey, mFragments.size(), finalised());
    }
}
