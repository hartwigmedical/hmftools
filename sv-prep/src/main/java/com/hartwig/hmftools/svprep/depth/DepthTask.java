package com.hartwig.hmftools.svprep.depth;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.svprep.reads.ReadGroup;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.variant.variantcontext.VariantContext;

public class DepthTask implements Callable
{
    private final List<VariantContext> mVariantsList;
    private int mBamRecords;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final Map<String,SamReadGroup> mReadGroups;

    public DepthTask(final String refGenomeFile, final String bamFile)
    {
        mVariantsList = Lists.newArrayList();
        mBamRecords = 0;

        mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(refGenomeFile)).open(new File(bamFile));
        mBamSlicer = new BamSlicer(0, true, true, true);
        mReadGroups = Maps.newHashMap();
    }

    public List<VariantContext> variants() { return mVariantsList; }

    @Override
    public Long call()
    {
        for(VariantContext variant : mVariantsList)
        {

        }



        return (long)0;
    }

    private void retrieveDepth(final VariantContext variant)
    {
        mReadGroups.clear();

        int variantPosition = variant.getStart();

        mBamSlicer.slice(mSamReader, Lists.newArrayList(new ChrBaseRegion(
                variant.getContig(), variantPosition, variantPosition)), this::processSamRecord);


        for(SamReadGroup readGroup : mReadGroups.values())
        {
            /*
            int readsPosMin = readGroup.reads().stream().mapToInt(x -> x.start()).min().orElse(0);
            int readsPosMax = readGroup.reads().stream().mapToInt(x -> x.end()).max().orElse(0);
            int baseStart = max(readsPosMin - mRegion.start(), 0);
            int baseEnd = min(readsPosMax - mRegion.start(), mBaseDepth.length - 1);
            for(int i = baseStart; i <= baseEnd; ++i)
            {
                ++mBaseDepth[i];
            }
            */
        }

    }

    private void processSamRecord(final SAMRecord record)
    {
        SamReadGroup readGroup = mReadGroups.get(record.getReadName());

        if(readGroup == null)
        {
            mReadGroups.put(record.getReadName(), new SamReadGroup(record));
            return;
        }

        readGroup.Reads.add(record);
    }

    private class SamReadGroup
    {
        public final List<SAMRecord> Reads;

        public SamReadGroup(final SAMRecord record)
        {
            Reads = Lists.newArrayList(record);
        }
    }

}
