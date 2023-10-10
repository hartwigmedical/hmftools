package com.hartwig.hmftools.sage.quality;

import static java.lang.Math.max;

import java.io.File;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.sage.ReferenceData;
import com.hartwig.hmftools.sage.SageConfig;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;

public class BamSampler
{
    private final SageConfig mConfig;
    private final ReferenceData mReferenceData;
    private final BamSlicer mSlicer;
    private int mReadCount;
    private int mMaxReadLength;

    private static final int MAX_READS = 1000;

    public BamSampler(final SageConfig config, final ReferenceData referenceData)
    {
        mConfig = config;
        mReferenceData = referenceData;
        mReadCount = 0;
        mMaxReadLength = 0;

        mSlicer = new BamSlicer(mConfig.MinMapQuality);
    }

    public void setBamCharacteristics(final String bamFile)
    {
        ChrBaseRegion sampleRegion = null;

        if(!mConfig.SpecificChrRegions.Regions.isEmpty())
        {
            sampleRegion = mConfig.SpecificChrRegions.Regions.get(0);
        }
        else if(!mReferenceData.PanelWithHotspots.isEmpty())
        {
            for(Map.Entry<Chromosome, List<BaseRegion>> entry : mReferenceData.PanelWithHotspots.entrySet())
            {
                BaseRegion region = entry.getValue().get(0);

                sampleRegion = new ChrBaseRegion(
                        mConfig.RefGenVersion.versionedChromosome(entry.getKey().toString()), region.start(), region.end());

                break;
            }
        }
        else
        {
            sampleRegion = new ChrBaseRegion(
                    mConfig.RefGenVersion.versionedChromosome(HumanChromosome._1.toString()), 100000, 200000);
        }


        SamReader samReader = SamReaderFactory.makeDefault()
                .referenceSource(new ReferenceSource(mReferenceData.RefGenome))
                .open(new File(bamFile));

        mSlicer.slice(samReader, sampleRegion, this::processRecord);

        if(mMaxReadLength > mConfig.getMaxObservedReadLength())
        {
            mConfig.setMaxObservedReadLength(mMaxReadLength);
        }
    }

    private void processRecord(final SAMRecord record)
    {
        ++mReadCount;

        mMaxReadLength = max(mMaxReadLength, record.getReadBases().length);

        if(mReadCount >= MAX_READS)
        {
            mSlicer.haltProcessing();
        }
    }
}
