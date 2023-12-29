package com.hartwig.hmftools.esvee.read;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.SvConstants.BAM_HEADER_SAMPLE_ID_TAG;
import static com.hartwig.hmftools.esvee.read.ReadCache.CACHE_QUERY_BUFFER;
import static com.hartwig.hmftools.esvee.read.ReadCache.MAX_CACHE_SIZE;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.esvee.SvConfig;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamReader
{
    private final SvConfig mConfig;

    private final List<SamReader> mSamReaders;
    private final BamSlicer mBamSlicer;

    private final ReadCache mReadCache;

    public BamReader(final SvConfig config)
    {
        mConfig = config;

        mSamReaders = Lists.newArrayList();

        for(int i = 0; i < config.SampleNames.size(); ++i)
        {
            String sampleId = config.SampleNames.get(i);
            String bamFile = config.BamFiles.get(i);

            SamReader samReader = SamReaderFactory.makeDefault()
                            .validationStringency(mConfig.BamStringency)
                            .referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(bamFile));

            samReader.getFileHeader().setAttribute(BAM_HEADER_SAMPLE_ID_TAG, sampleId);

            mSamReaders.add(samReader);
        }

        mBamSlicer = createBamSlicer();

        mReadCache = new ReadCache(MAX_CACHE_SIZE, new ReadCache.CacheLoader()
        {
            @Override
            public ReadCache.CachedReadKey determineCacheRange(final String chromosome, final int positionStart, final int positionEnd)
            {
                final int adjustedStart = Math.max(1, positionStart - CACHE_QUERY_BUFFER);
                final int adjustedEnd = positionEnd + CACHE_QUERY_BUFFER;

                return new ReadCache.CachedReadKey(chromosome, adjustedStart, adjustedEnd);
            }

            @Override
            public ReadCache.CachedReadValue load(final ReadCache.CachedReadKey key)
            {
                List<Read> alignments = doSliceBam(key.Chromosome, key.PositionStart, key.PositionEnd);
                return new ReadCache.CachedReadValue(alignments);
            }
        });
    }

    public void sliceBam(final String chromosome, int positionStart, int positionEnd, final Consumer<SAMRecord> consumer)
    {
        int bamPosStart = max(positionStart, 1);
        int bamPosEnd = min(positionEnd, mConfig.RefGenomeCoords.length(chromosome));

        if(bamPosStart > bamPosEnd)
            return;

        for(SamReader reader : mSamReaders)
        {
            mBamSlicer.slice(reader, new ChrBaseRegion(chromosome, positionStart, positionEnd), consumer);
        }
    }

    public List<Read> sliceBam(final String chromosome, int positionStart, int positionEnd)
    {
        int bamPosStart = max(positionStart, 1);
        int bamPosEnd = min(positionEnd, mConfig.RefGenomeCoords.length(chromosome));

        if(bamPosStart > bamPosEnd)
            return Collections.emptyList();

        try
        {
            return mReadCache.read(chromosome, bamPosStart, bamPosEnd).collect(Collectors.toList());
        }
        finally
        {
            /* FIXME: decide how to log and measure effectiveness
            int hits = mReadCache.CacheHits.get();
            int misses = mReadCache.CacheMisses.get();

            if((hits + misses) % 10_000 == 0)
            {
                SV_LOGGER.info("read cache: cache hits({}) misses({}) hitRate({}%)",
                        hits, misses, String.format("%.2f", (100.0f * hits) / (hits + misses)));
            }
            */
        }
    }

    private List<Read> doSliceBam(final String chromosome, int bamPosStart, int bamPosEnd)
    {
        ChrBaseRegion sliceRegion = new ChrBaseRegion(chromosome, bamPosStart, bamPosEnd);

        List<Read> reads = Lists.newArrayList();

        for(SamReader reader : mSamReaders)
        {
            // FIXME: consumer or accumulator
            final Consumer<SAMRecord> consumer = record ->
            {
                reads.add(createAndCheckRecord(record));
            };

            mBamSlicer.slice(reader, sliceRegion, consumer);
        }

        return reads;
    }

    private Read createAndCheckRecord(final SAMRecord record)
    {
        // mReadRescue::rescueRead) // CHECK: not required

        // mNormaliser::normalise() on each record

        return new Read(record);
    }

    private BamSlicer createBamSlicer()
    {
        // as per SvPrep
        BamSlicer bamSlicer = new BamSlicer(0, false, true, false);
        bamSlicer.setKeepUnmapped();
        bamSlicer.setKeepHardClippedSecondaries();
        return bamSlicer;
    }

}
