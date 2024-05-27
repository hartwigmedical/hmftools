package com.hartwig.hmftools.bamtools.compare;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.bamtools.common.HashBamWriter;

import org.apache.commons.lang3.tuple.Pair;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class UnmatchedReadHandler
{
    private static final int NUM_HASH_BAMS = 512;

    private final HashBamWriter mOrigBamHashBamWriter;
    private final HashBamWriter mNewBamHashBamWriter;

    public UnmatchedReadHandler(final CompareConfig config)
    {
        try(SamReader samReader = CompareUtils.makeSamReaderFactory(config).open(new File(config.OrigBamFile)))
        {
            SAMFileHeader header = samReader.getFileHeader();
            // create hash bams
            mOrigBamHashBamWriter = new HashBamWriter(header, "bamcomp_orig_hashbams_", NUM_HASH_BAMS);
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }

        try(SamReader samReader = CompareUtils.makeSamReaderFactory(config).open(new File(config.NewBamFile)))
        {
            SAMFileHeader header = samReader.getFileHeader();
            // create hash bams
            mNewBamHashBamWriter = new HashBamWriter(header, "bamcomp_new_hashbams_", NUM_HASH_BAMS);
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }

    void closeHashBamWriters()
    {
        mOrigBamHashBamWriter.close();
        mNewBamHashBamWriter.close();
    }

    public void handleOrigBamReads(Collection<SAMRecord> reads)
    {
        BT_LOGGER.trace("writing {} orig bam reads to hash bams", reads.size());
        writeReadsToHashBams(mOrigBamHashBamWriter, reads.iterator());
    }
    public void handleNewBamReads(Collection<SAMRecord> reads)
    {
        BT_LOGGER.trace("writing {} new bam reads to hash bams", reads.size());
        writeReadsToHashBams(mNewBamHashBamWriter, reads.iterator());
    }

    public long handleOrigBamReads(Iterator<SAMRecord> itr)
    {
        return writeReadsToHashBams(mOrigBamHashBamWriter, itr);
    }

    public long handleNewBamReads(Iterator<SAMRecord> itr)
    {
        return writeReadsToHashBams(mNewBamHashBamWriter, itr);
    }

    private static long writeReadsToHashBams(HashBamWriter hashBamWriter, Iterator<SAMRecord> itr)
    {
        long numReads = 0;
        while(itr.hasNext())
        {
            hashBamWriter.writeToHashBam(itr.next());
            ++numReads;
        }
        return numReads;
    }

    // get the key -> hashbam pairs
    public Map<Integer, Pair<File, File>> getHashBamPairs()
    {
        // get all the keys, and double check they are correct
        Map<Integer, File> origBamHashBams = mOrigBamHashBamWriter.getHashBams();
        Map<Integer, File> newBamHashBams = mNewBamHashBamWriter.getHashBams();
        Set<Integer> hashBamKeys = new HashSet<>(origBamHashBams.keySet());
        hashBamKeys.addAll(newBamHashBams.keySet());

        Map<Integer, Pair<File, File>> hashBamPairs = new HashMap<>();

        for(Integer key : hashBamKeys)
        {
            File origHashBam = origBamHashBams.get(key);
            File newHashBam = newBamHashBams.get(key);

            // should not be null
            if(origHashBam == null || newHashBam == null)
            {
                throw new RuntimeException(String.format("missing hashbams with key: %d", key));
            }

            hashBamPairs.put(key, Pair.of(origHashBam, newHashBam));
        }
        return hashBamPairs;
    }
}
