package com.hartwig.hmftools.redux.common;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;

public class PartitionResults
{
    private List<Fragment> mResolvedFragments;
    private List<DuplicateGroup> mDuplicateGroups;
    private FragmentStatus mFragmentStatus; // set for a single non-primary fragment

    public PartitionResults()
    {
        mResolvedFragments = null;
        mDuplicateGroups = null;
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

    public List<DuplicateGroup> umiGroups() { return mDuplicateGroups; }

    public void addUmiGroup(final DuplicateGroup umiGroup)
    {
        if(mDuplicateGroups == null)
            mDuplicateGroups = Lists.newArrayList(umiGroup);
        else
            mDuplicateGroups.add(umiGroup);
    }

    public String toString()
    {
        return format("fragStatus(%s) resolved(%d) umiGroups(%d)",
                mFragmentStatus != null ? mFragmentStatus : "unset",
                mResolvedFragments != null ? mResolvedFragments.size() : "unset",
                mDuplicateGroups != null ? mDuplicateGroups.size() : "unset");
    }

}
