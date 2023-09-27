package com.hartwig.hmftools.sieve.annotate;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.sieve.annotate.AnnotateConfig.MD_LOGGER;
import static com.hartwig.hmftools.sieve.annotate.Util.getNM;
import static com.hartwig.hmftools.sieve.annotate.Util.getSoftClipCount;
import static com.hartwig.hmftools.sieve.annotate.Util.isDiscordantPair;

import htsjdk.samtools.SAMRecord;

public class HighDepthCounts
{
    // Bucket MapQs into 0-9, 10-19, ..., 40-49, >= 50.
    private static final int MAPQ_BUCKET_COUNT = 6;
    private static final int NM_BUCKET_COUNT = 11;
    // Bucket SOFT_CLIP counts into buckets of 1-3, 4-10, and >=11.
    private static final int SOFT_CLIP_BUCKET_COUNT = 3;
    public static final String TSV_HEADER = getTSVHeader();

    private final long[] mConcordantMapQBuckets;
    private final long[] mConcordantNMBuckets;
    private final long[] mConcordantSoftClipBuckets;

    private long mSupplementaryCount;
    private long mPrimaryReadCount;
    private long mPrimaryDiscordantCount;
    private long mPrimaryConcordantCount;
    private long mPrimaryConcordantShortReads;

    public HighDepthCounts()
    {
        mSupplementaryCount = 0;
        mPrimaryReadCount = 0;
        mPrimaryDiscordantCount = 0;
        mPrimaryConcordantCount = 0;
        mPrimaryConcordantShortReads = 0;

        mConcordantMapQBuckets = new long[MAPQ_BUCKET_COUNT];
        for(int i = 0; i < mConcordantMapQBuckets.length; i++)
        {
            mConcordantMapQBuckets[i] = 0;
        }

        mConcordantNMBuckets = new long[NM_BUCKET_COUNT];
        for(int i = 0; i < mConcordantNMBuckets.length; i++)
        {
            mConcordantNMBuckets[i] = 0;
        }

        mConcordantSoftClipBuckets = new long[SOFT_CLIP_BUCKET_COUNT];
        for(int i = 0; i < mConcordantSoftClipBuckets.length; i++)
        {
            mConcordantSoftClipBuckets[i] = 0;
        }
    }

    private static String getTSVHeader()
    {
        final StringBuilder sb = new StringBuilder();
        sb.append("SupplementaryCount");
        sb.append('\t');
        sb.append("PrimaryReadCount");
        sb.append('\t');
        sb.append("PrimaryDiscordantCount");
        sb.append('\t');
        sb.append("PrimaryConcordantCount");
        sb.append('\t');
        sb.append("PrimaryConcordantShortReads");

        for(int i = 0; i < MAPQ_BUCKET_COUNT; i++)
        {
            final int mapQStart = 10 * i;
            final int mapQEnd = mapQStart + 9;
            if(i < MAPQ_BUCKET_COUNT - 1)
            {
                sb.append("\tConcordant MAPQ " + mapQStart + '-' + mapQEnd);
            }
            else
            {
                sb.append("\tConcordant MAPQ >=" + mapQStart);
            }
        }

        for(int i = 0; i < NM_BUCKET_COUNT; i++)
        {
            if(i < NM_BUCKET_COUNT - 1)
            {
                sb.append("\tConcordant NM " + i);
            }
            else
            {
                sb.append("\tConcordant NM >=" + i);
            }
        }

        // Bucket SOFT_CLIP counts into buckets of 1-3, 4-10, and >=11.
        sb.append("\tConcordant SOFT_CLIP 1-3");
        sb.append("\tConcordant SOFT_CLIP 4-10");
        sb.append("\tConcordant SOFT_CLIP >=11");

        return sb.toString();
    }

    private static int getMapQBucketIndex(int mapQ)
    {
        return Math.min(mapQ / 10, MAPQ_BUCKET_COUNT - 1);
    }

    private static int getNMBucketIndex(int nm)
    {
        return Math.min(nm, NM_BUCKET_COUNT - 1);
    }

    private static int getSoftClipBucketIndex(int softClipSize)
    {
        if(softClipSize <= 0)
        {
            MD_LOGGER.error("Cannot get soft clip bucket index when there are no soft clips.");
            System.exit(1);
        }

        // Bucket SOFT_CLIP counts into buckets of 1-3, 4-10, and >=11.
        if(softClipSize <= 3)
        {
            return 0;
        }

        if(softClipSize <= 10)
        {
            return 1;
        }

        return 2;
    }

    public void matchedRead(final SAMRecord read)
    {
        if(read.getSupplementaryAlignmentFlag())
        {
            mSupplementaryCount++;
            return;
        }

        mPrimaryReadCount++;

        if(isDiscordantPair(read))
        {
            mPrimaryDiscordantCount++;
            return;
        }

        mPrimaryConcordantCount++;
        // TODO(m_cooper): Make 150 configurable?
        if(abs(read.getInferredInsertSize()) < 150)
        {
            mPrimaryConcordantShortReads++;
        }

        mConcordantMapQBuckets[getMapQBucketIndex(read.getMappingQuality())]++;

        Integer editDistance = getNM(read);
        if(editDistance == null)
        {
            MD_LOGGER.error("Concordant primary read {} does not contain NM attribute.", read);
            System.exit(1);
        }

        mConcordantNMBuckets[getNMBucketIndex(editDistance)]++;

        int softClipCount = getSoftClipCount(read);
        if(softClipCount > 0)
        {
            mConcordantSoftClipBuckets[getSoftClipBucketIndex(softClipCount)]++;
        }
    }

    public String getTSVFragment()
    {
        StringBuilder sb = new StringBuilder();
        sb.append(String.valueOf(mSupplementaryCount) + '\t' + mPrimaryReadCount + '\t' + mPrimaryDiscordantCount + '\t'
                + mPrimaryConcordantCount + '\t' + mPrimaryConcordantShortReads);

        for(int i = 0; i < mConcordantMapQBuckets.length; i++)
        {
            sb.append('\t');
            sb.append(mConcordantMapQBuckets[i]);
        }

        for(int i = 0; i < mConcordantNMBuckets.length; i++)
        {
            sb.append('\t');
            sb.append(mConcordantNMBuckets[i]);
        }

        for(int i = 0; i < mConcordantSoftClipBuckets.length; i++)
        {
            sb.append('\t');
            sb.append(mConcordantSoftClipBuckets[i]);
        }

        return sb.toString();
    }
}
