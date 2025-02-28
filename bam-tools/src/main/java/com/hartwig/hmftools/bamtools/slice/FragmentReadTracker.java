package com.hartwig.hmftools.bamtools.slice;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.bam.SamRecordUtils;
import htsjdk.samtools.SAMRecord;

import org.apache.logging.log4j.Level;

import java.util.ArrayList;
import java.util.List;

import javax.annotation.Nullable;

// The FragmentReadTracker class is responsible for tracking the status of reads
// in a sequencing fragment. It processes single-end or paired-end reads along
// with their supplementary alignments, ensuring all reads in a fragment are accounted for.
// It maintains a list of pending and processed reads and provides utilities for:
// - Adding and processing reads.
// - Determining if a fragment is complete.
// - Finding genomic positions of pending reads for further processing.
//
// Assumptions:
// - Mate position: If the read is paired, mate position is taken from the SAMRecord, except in the case
//   where current alignment is supplementary and the mate is unmapped, in this case the mate alignment
//   is derived from the first alignment in the SA (supplementary alignment) tag.
// - SA tag: SA tag is populated. For supplementary reads, the first entry in SA tag is the primary read.
//
public class FragmentReadTracker
{
    public record ReadKey(boolean isSupplementary, boolean firstInPair, boolean isNegativeStrand,
                          boolean isSecondary, int alignmentStart,
                          String referenceName, @Nullable String cigar)
    {
        public ReadKey {
            if(isSupplementary && cigar == null)
            {
                throw new IllegalArgumentException("cigar string missing for supplementary read");
            }
            if(!isSupplementary && cigar != null)
            {
                throw new IllegalArgumentException("cigar string present for primary read");
            }
        }

        @Override
        public String toString()
        {
            return String.format("supplementary?=%b, firstInPair?=%b, isSecondary?=%b, aligned to %s:%d(%c), cigar=%s",
                    isSupplementary,
                    firstInPair,
                    isSecondary,
                    referenceName,
                    alignmentStart,
                    isNegativeStrand ? '-' : '+',
                    cigar);
        }

        static ReadKey fromRead(SAMRecord read)
        {
            // NOTE: we only need cigar string for supplementary reads
            // also use intern String to save memory as most cigar strings are 151M
            String cigarString = read.getSupplementaryAlignmentFlag() ? read.getCigarString().intern() : null;
            return new ReadKey(
                    read.getSupplementaryAlignmentFlag(),
                    SamRecordUtils.firstInPair(read),
                    read.getReadNegativeStrandFlag(),
                    read.isSecondaryAlignment(),
                    read.getAlignmentStart(),
                    read.getReferenceName(),
                    cigarString);
        }
    }

    private final String mName;
    private final boolean mIsReadPaired;
    private final List<ReadKey> mPendingReads = new ArrayList<>();
    private final List<ReadKey> mProcessedReads = new ArrayList<>();

    public FragmentReadTracker(String name, boolean isReadPaired)
    {
        mName = name;
        mIsReadPaired = isReadPaired;
    }

    public String getName()
    {
        return mName;
    }

    public List<ReadKey> getPendingReads()
    {
        return mPendingReads;
    }

    public boolean isComplete(Level logLevel)
    {
        if(!mPendingReads.isEmpty())
        {
            if(BT_LOGGER.isEnabled(logLevel))
            {
                for(ReadKey readKey : mPendingReads)
                {
                    BT_LOGGER.log(logLevel, "read id({}) pending read({}) ", mName, readKey);
                }
            }
            return false;
        }
        if(mIsReadPaired)
        {
            // paired reads then we need minimum of two non supplementary
            if(mProcessedReads.size() < 2 || mProcessedReads.stream().filter(x -> !x.isSupplementary).count() < 2)
            {
                return false;
            }
        }
        else
        {
            if(mProcessedReads.isEmpty() || mProcessedReads.stream().allMatch(x -> x.isSupplementary))
            {
                return false;
            }
        }

        // if we got all the primary reads, that means we got all the supplementary information as well
        // and then if pending reads is not empty then it is complete
        return true;
    }

    public boolean invariant()
    {
        if(mIsReadPaired && mPendingReads.isEmpty() && mProcessedReads.stream().filter(x -> !x.isSupplementary).count() != 2)
        {
            // if read is paired then we expect to have two processed reads
            return false;
        }
        if(!mIsReadPaired && mPendingReads.isEmpty() && mProcessedReads.stream().filter(x -> !x.isSupplementary).count() != 1)
        {
            // if read is not paired then we expect to have one processed reads
            return false;
        }
        return true;
    }

    public List<BasePosition> findPendingReadPositions()
    {
        List<BasePosition> baseRegions = new ArrayList<>();
        assert invariant();

        for(ReadKey readKey : mPendingReads)
        {
            baseRegions.add(new BasePosition(readKey.referenceName, readKey.alignmentStart));
        }
        return baseRegions;
    }

    public boolean processRead(SAMRecord read)
    {
        ReadKey readKey = ReadKey.fromRead(read);

        if(mProcessedReads.stream().anyMatch(x -> x.equals(readKey)))
        {
            // we have already seen this
            return false;
        }

        if(read.isSecondaryAlignment())
        {
            // secondary alignments are added the processed reads to make sure they are recorded, but
            // they are not used to find other pending alignments
            mProcessedReads.add(readKey);
            return true;
        }

        // add pending mate + supplementary alignments
        // we only need to do this for the first alignment of the read
        if(mProcessedReads.stream().noneMatch(o -> o.firstInPair == SamRecordUtils.firstInPair(read)))
        {
            // we want to find all the supplementary alignments and add to pending
            List<SupplementaryReadData> supplementaryAlignments = SupplementaryReadData.extractAlignments(read);
            if(supplementaryAlignments != null)
            {
                for(int i = 0; i < supplementaryAlignments.size(); ++i)
                {
                    SupplementaryReadData sa = supplementaryAlignments.get(i);

                    // if this read is supplementary, then first read in SA tag is primary
                    boolean isSupplementary = !read.getSupplementaryAlignmentFlag() || i > 0;
                    ReadKey saReadKey = new ReadKey(isSupplementary,
                            SamRecordUtils.firstInPair(read),
                            sa.Strand == SupplementaryReadData.SUPP_NEG_STRAND,
                            false,
                            sa.Position,
                            sa.Chromosome,
                            isSupplementary ? sa.Cigar.intern() : null);

                    if(mPendingReads.stream().noneMatch(saReadKey::equals))
                    {
                        // no need to check if processedReads contains this saReadKey as the outer
                        // if case protects against it
                        mPendingReads.add(saReadKey);
                        BT_LOGGER.trace("{} added SA read key({}) to pending", read, saReadKey);
                    }
                }
            }

            // add the mate read key
            // NOTE: if current alignment is supplementary, the mate is still the primary alignment
            // only problem is if this alignment is supplementary and mate is unmapped, then the mate alignment
            // will show this supplementary alignment instead, which is not what we want. We want to extract
            // the primary alignment instead from the SA tag
            if(read.getReadPairedFlag())
            {
                ReadKey mateReadKey = null;
                if(!read.getSupplementaryAlignmentFlag() || !read.getMateUnmappedFlag())
                {
                    mateReadKey = new ReadKey(false,
                            !SamRecordUtils.firstInPair(read),
                            read.getMateNegativeStrandFlag(),
                            false,
                            read.getMateAlignmentStart(),
                            read.getMateReferenceName(),
                            null); // do not populate cigar for primary alignment
                }
                else if(supplementaryAlignments != null && !supplementaryAlignments.isEmpty())
                {
                    // if this alignment is supplementary and mate is unmapped, we cannot simply take the
                    // mate alignment from this SAMRecord, it will show this supplementary alignment which is incorrect.
                    // Instead we have to use the first alignment in the SA tag
                    SupplementaryReadData sa = supplementaryAlignments.get(0);
                    mateReadKey = new ReadKey(false,
                            !SamRecordUtils.firstInPair(read),
                            sa.Strand == SupplementaryReadData.SUPP_NEG_STRAND,
                            false,
                            sa.Position,
                            sa.Chromosome,
                            null); // do not populate cigar for primary alignment
                }
                if(mateReadKey != null && mProcessedReads.stream().noneMatch(mateReadKey::equals) &&
                        mPendingReads.stream().noneMatch(mateReadKey::equals))
                {
                    mPendingReads.add(mateReadKey);
                    BT_LOGGER.trace("{} added mate read key({}) to pending", read, mateReadKey);
                }
            }
        }
        mProcessedReads.add(readKey);
        mPendingReads.removeIf(x -> x.equals(readKey));

        return true;
    }
}
