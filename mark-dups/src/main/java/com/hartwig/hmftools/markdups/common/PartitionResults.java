package com.hartwig.hmftools.markdups.common;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.markdups.umi.UmiGroup;

public class PartitionResults
{
    private List<Fragment> mResolvedFragments;
    private List<UmiGroup> mUmiGroups;
    private FragmentStatus mFragmentStatus; // set for a single non-primary fragment

    public PartitionResults()
    {
        mResolvedFragments = null;
        mUmiGroups = null;
        mFragmentStatus = null;
    }

    public FragmentStatus fragmentStatus() { return mFragmentStatus; }
    public void setFragmentStatus(final FragmentStatus status) { mFragmentStatus = status; }

    public List<Fragment> resolvedFragments() { return mResolvedFragments; }

    public void addResolvedFragments(final List<Fragment> fragments)
    {
        if(mResolvedFragments == null)
            mResolvedFragments = Lists.newArrayList(fragments);
        else
            mResolvedFragments.addAll(fragments);
    }

    public void addResolvedFragment(final Fragment fragment)
    {
        if(mResolvedFragments == null)
            mResolvedFragments = Lists.newArrayList(fragment);
        else
            mResolvedFragments.add(fragment);
    }

    public List<UmiGroup> umiGroups() { return mUmiGroups; }

    public void addUmiGroup(final UmiGroup umiGroup)
    {
        if(mUmiGroups == null)
            mUmiGroups = Lists.newArrayList(umiGroup);
        else
            mUmiGroups.add(umiGroup);
    }

    public String toString() { return format("resolved(%d) umiGroups(%d)", mResolvedFragments.size(), mUmiGroups.size()); }

}
