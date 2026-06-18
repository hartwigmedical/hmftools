package com.hartwig.hmftools.tars.liftback;

import java.util.HashMap;
import java.util.Map;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;

// Maps a lifted record to a genomic bucket so the shards can be sorted per bucket and then plain-concatenated
// (cat) in bucket order to produce a coordinate-sorted BAM -- buckets are disjoint, ordered ranges of the
// absolute genome offset, so cat of the sorted pieces is globally sorted. Routing is by reference NAME (the
// record may still carry the input header at write time, so its reference index is not reliable). Unmapped
// records (reference "*") go to the final bucket, matching coordinate sort's placement of them at the end.
public class GenomicBucketIndex
{
    private final Map<String,Long> mContigOffset;
    private final long mBinWidth;
    private final int mMappedBuckets;

    public GenomicBucketIndex(final SAMFileHeader header, final int mappedBuckets)
    {
        mMappedBuckets = Math.max(1, mappedBuckets);
        mContigOffset = new HashMap<>();

        long offset = 0;
        for(final SAMSequenceRecord seq : header.getSequenceDictionary().getSequences())
        {
            mContigOffset.put(seq.getSequenceName(), offset);
            offset += seq.getSequenceLength();
        }

        final long totalLength = Math.max(1, offset);
        mBinWidth = (totalLength + mMappedBuckets - 1) / mMappedBuckets;
    }

    // total bucket count: mapped bins plus one trailing bucket for unmapped records.
    public int bucketCount() { return mMappedBuckets + 1; }

    public int bucketOf(final SAMRecord record)
    {
        final String chromosome = record.getReferenceName();
        if(chromosome == null || chromosome.equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME))
            return mMappedBuckets;

        final Long offset = mContigOffset.get(chromosome);
        if(offset == null)
            throw new IllegalStateException("lifted record references contig absent from output header: " + chromosome);

        final long absolutePosition = offset + record.getAlignmentStart();
        final int bucket = (int)(absolutePosition / mBinWidth);
        return Math.min(bucket, mMappedBuckets - 1);
    }
}
