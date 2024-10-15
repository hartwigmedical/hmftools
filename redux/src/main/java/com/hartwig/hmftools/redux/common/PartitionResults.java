package com.hartwig.hmftools.redux.common;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;

import htsjdk.samtools.SAMRecord;

public class PartitionResults
{
    private List<Fragment> mResolvedFragments;
    private List<DuplicateGroup> mDuplicateGroups;
    private FragmentStatus mFragmentStatus; // set for a single non-primary fragment
    private List<SAMRecord> mSupplementaries;

    public PartitionResults()
    {
        mResolvedFragments = null;
        mDuplicateGroups = null;
        mFragmentStatus = null;
        mSupplementaries = null;
    }

    public PartitionResults(final SAMRecord supplementaryRead)
    {
        mResolvedFragments = null;
        mDuplicateGroups = null;
        mFragmentStatus = null;
        mSupplementaries = List.of(supplementaryRead);
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

    public List<SAMRecord> supplementaries() { return mSupplementaries; }

    public void addSupplementary(final SAMRecord read)
    {
        if(mSupplementaries == null)
            mSupplementaries = Lists.newArrayList(read);
        else
            mSupplementaries.add(read);
    }

    public String toString()
    {
        return format("fragStatus(%s) resolved(%d) umiGroups(%d) supplementaries(%d)",
                mFragmentStatus != null ? mFragmentStatus : "unset",
                mResolvedFragments != null ? mResolvedFragments.size() : 0,
                mDuplicateGroups != null ? mDuplicateGroups.size() : 0,
                mSupplementaries != null ? mSupplementaries.size() : 0);
    }

}
