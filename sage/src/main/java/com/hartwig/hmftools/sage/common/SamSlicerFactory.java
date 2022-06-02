package com.hartwig.hmftools.sage.common;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.SageConfig;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class SamSlicerFactory
{
    private final Map<String,SamSlicerInterface> mSamSlicers;

    private final Map<String,SamReader> mBamReaders;

    public SamSlicerFactory()
    {
        mSamSlicers = Maps.newHashMap();
        mBamReaders = Maps.newHashMap();
    }

    public SamSlicerInterface getSamSlicer(final String sampleId, final List<ChrBaseRegion> regions)
    {
        if(!mBamReaders.isEmpty())
        {
            SamReader bamReader = mBamReaders.get(sampleId);
            return new SamSlicer(bamReader, 0, regions);
        }

        return mSamSlicers.get(sampleId);
    }

    public void buildBamReaders(final SageConfig config, final IndexedFastaSequenceFile refGenome)
    {
        List<String> allSamples = Lists.newArrayList(config.TumorIds);
        allSamples.addAll(config.ReferenceIds);

        List<String> allBams = Lists.newArrayList(config.TumorBams);
        allBams.addAll(config.ReferenceBams);

        for(int i = 0; i < allSamples.size(); i++)
        {
            final String sample = allSamples.get(i);
            final String bamFile = allBams.get(i);

            SamReader bamReader = SamReaderFactory.makeDefault()
                    .validationStringency(config.Stringency)
                    .referenceSource(new ReferenceSource(refGenome))
                    .open(new File(bamFile));

            mBamReaders.put(sample, bamReader);
        }
    }

    public void close()
    {
        try
        {
            for(SamReader bamReader : mBamReaders.values())
            {
                bamReader.close();
            }
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to close bam file: {}", e.toString());
        }
    }
}
