package com.hartwig.hmftools.bamtools.slice;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.bam.SamRecordUtils;
import htsjdk.samtools.SAMRecord;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.ArrayList;
import java.util.List;

import javax.annotation.Nullable;

public class FragmentReadTracker
{
    private static final Logger logger = LogManager.getLogger(FragmentReadTracker.class);

    public record ReadKey(boolean secondaryOrSupplementary, boolean firstInPair, boolean isNegativeStrand,
                          int alignmentStart,
                          String referenceName, @Nullable String cigar)
    {
        public ReadKey {
            if(secondaryOrSupplementary && cigar == null)
            {
                throw new IllegalArgumentException("cigar string missing for supplementary read");
            }
            if(!secondaryOrSupplementary && cigar != null)
            {
                throw new IllegalArgumentException("cigar string present for primary read");
            }
        }

        @Override
        public String toString()
        {
            return "2ndaryOrSuppl?=" + secondaryOrSupplementary + ", firstInPair?=" + firstInPair +
                    ", aligned to " + referenceName + ':' + alignmentStart + '(' +
                    (isNegativeStrand ? '-' : '+') + "), cigar=" + cigar;
        }

        static ReadKey fromRead(SAMRecord read)
        {
            // NOTE: we only need cigar string for supplementary reads
            // also use intern String to save memory as most cigar strings are 151M
            String cigarString = read.isSecondaryOrSupplementary() ? read.getCigarString().intern() : null;
            return new ReadKey(
                    read.isSecondaryOrSupplementary(),
                    SamRecordUtils.firstInPair(read),
                    read.getReadNegativeStrandFlag(),
                    read.getAlignmentStart(),
                    read.getReferenceName(),
                    cigarString);
        }
    }

    private final String name;
    private final boolean isReadPaired;
    private final List<ReadKey> pendingReads = new ArrayList<>();
    private final List<ReadKey> processedReads = new ArrayList<>();

    public FragmentReadTracker(String name, boolean isReadPaired)
    {
        this.name = name;
        this.isReadPaired = isReadPaired;
    }

    public String getName()
    {
        return name;
    }

    public List<ReadKey> getPendingReads()
    {
        return pendingReads;
    }

    public boolean isComplete(Level logLevel)
    {
        if(!pendingReads.isEmpty())
        {
            if(logger.isEnabled(logLevel))
            {
                for(ReadKey readKey : pendingReads)
                {
                    logger.log(logLevel, "read id({}) pending read({}) ", name, readKey);
                }
            }
            return false;
        }
        if(isReadPaired)
        {
            // paired reads then we need minimum of two non supplementary
            if(processedReads.size() < 2 || processedReads.stream().filter(x -> !x.secondaryOrSupplementary).count() < 2)
            {
                return false;
            }
        }
        else
        {
            if(processedReads.isEmpty() || processedReads.stream().allMatch(x -> x.secondaryOrSupplementary))
                return false;
        }

        // if we got all the primary reads, that means we got all the supplementary information as well
        // and then if pending reads is not empty then it is complete
        return true;
    }

    @Override
    public String toString()
    {
        return String.format("%s processed reads(%d), pending reads(%d) complete(%s)",
                name, processedReads.size(), pendingReads.size(), isComplete(null));
    }

    public boolean invariant()
    {
        if(isReadPaired && pendingReads.isEmpty() && processedReads.stream().filter(x -> !x.secondaryOrSupplementary).count() != 2)
        {
            // if read is paired then we expect to have two processed reads
            return false;
        }
        if(!isReadPaired && pendingReads.isEmpty() && processedReads.stream().filter(x -> !x.secondaryOrSupplementary).count() != 1)
        {
            // if read is not paired then we expect to have one processed reads
            return false;
        }
        return true;
    }

    public List<BasePosition> findPendingReadBaseRegions()
    {
        List<BasePosition> baseRegions = new ArrayList<>();
        assert invariant();

        for(ReadKey readKey : pendingReads)
        {
            baseRegions.add(new BasePosition(readKey.referenceName, readKey.alignmentStart));
        }
        return baseRegions;
    }

    public boolean processRead(SAMRecord read)
    {
        ReadKey readKey = ReadKey.fromRead(read);

        if(processedReads.stream().anyMatch(x -> x.equals(readKey)))
        {
            // we have already seen this
            return false;
        }

        // add pending mate + supplementary alignments
        // we only need to do this for the first alignment of the read
        if(processedReads.stream().noneMatch(o -> o.firstInPair == SamRecordUtils.firstInPair(read)))
        {
            // we want to find all the supplementary alignments and add to pending
            List<SupplementaryReadData> supplementaryAlignments = SupplementaryReadData.extractAlignments(read);
            if(supplementaryAlignments != null)
            {
                for(int i = 0; i < supplementaryAlignments.size(); ++i)
                {
                    SupplementaryReadData sa = supplementaryAlignments.get(i);

                    // if this read is supplementary, then first read in SA tag is primary
                    boolean isSupplementary = !read.isSecondaryOrSupplementary() || i > 0;
                    ReadKey saReadKey = new ReadKey(isSupplementary,
                            SamRecordUtils.firstInPair(read),
                            sa.Strand == SupplementaryReadData.SUPP_NEG_STRAND,
                            sa.Position,
                            sa.Chromosome,
                            isSupplementary ? sa.Cigar.intern() : null);

                    if(pendingReads.stream().noneMatch(saReadKey::equals))
                    {
                        pendingReads.add(saReadKey);
                        // logger.trace("{} Missing supplementary read: aligned to {}:{}", read, sa.Chromosome, sa.Position);
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
                if(!read.isSecondaryOrSupplementary() || !read.getMateUnmappedFlag())
                {
                    mateReadKey = new ReadKey(false,
                            !SamRecordUtils.firstInPair(read),
                            read.getMateNegativeStrandFlag(),
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
                            sa.Position,
                            sa.Chromosome,
                            null); // do not populate cigar for primary alignment
                }
                if(mateReadKey != null && processedReads.stream().noneMatch(mateReadKey::equals))
                {
                    pendingReads.add(mateReadKey);
                    logger.trace("{} added mate read key for pending: aligned to {}:{}",
                            read,
                            read.getMateReferenceName(),
                            read.getMateAlignmentStart());
                }
            }
        }
        processedReads.add(readKey);
        pendingReads.removeIf(x -> x.equals(readKey));

        return true;
    }
}
