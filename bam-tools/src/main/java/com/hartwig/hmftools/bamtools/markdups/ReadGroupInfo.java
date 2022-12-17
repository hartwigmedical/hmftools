package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.markdups.GroupCombiner.formChromosomePartition;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bamtools.ReadGroup;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

import htsjdk.samtools.SAMRecord;

public class ReadGroupInfo
{
    public final String ReadId;
    public final DuplicateStatus Status;
    public final boolean IsComplete;

    public final BaseRegion CurrentRange;
    public final BaseRegion ExpectedRange;
    public final List<String> ChrPartitions;

    public ReadGroupInfo(final ReadGroup readGroup, DuplicateStatus dupStatus, final BaseRegion currentPartition)
    {
        ReadId = readGroup.id();
        Status = dupStatus;

        List<String> chrPartitions = Lists.newArrayList();
        String chromosome = readGroup.reads().get(0).getContig();
        int partitionSize = currentPartition.baseLength();

        int expectedNonSuppCount = 1;
        int expectedSuppCount = 0;
        int nonSuppCount = 0;
        int suppCount = 0;

        BaseRegion currentRange = null;
        BaseRegion expectedRange = null;

        for(int i = 0; i < readGroup.reads().size(); ++i)
        {
            SAMRecord read = readGroup.reads().get(i);
            int readPosStart = read.getAlignmentStart();

            if(i == 0)
            {
                currentRange = new BaseRegion(readPosStart, readPosStart);
                expectedRange = new BaseRegion(readPosStart, readPosStart);
            }
            else
            {
                currentRange.setStart(min(currentRange.start(), readPosStart));
                currentRange.setEnd(max(currentRange.end(), readPosStart));
                expectedRange.setStart(min(expectedRange.start(), currentRange.start()));
                expectedRange.setEnd(max(expectedRange.end(), currentRange.end()));
            }

            if(read.getReadPairedFlag() && !read.getMateUnmappedFlag())
            {
                expectedNonSuppCount = 2;

                if(read.getMateReferenceName().equals(chromosome)) //  && currentPartition.containsPosition(read.getMateAlignmentStart())
                {
                    expectedRange.setStart(min(expectedRange.start(), read.getMateAlignmentStart()));
                    expectedRange.setEnd(max(expectedRange.end(), read.getMateAlignmentStart()));
                }
                else
                {
                    // a remote partition - do these really need to be recorded if the mate will also come to the same classification?
                    chrPartitions.add(formChromosomePartition(read.getMateReferenceName(), read.getMateAlignmentStart(), partitionSize));
                }
            }

            if(read.getSupplementaryAlignmentFlag())
            {
                ++suppCount;
            }
            else
            {
                ++nonSuppCount;
            }

            if(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
            {
                if(!read.getSupplementaryAlignmentFlag())
                {
                    ++expectedSuppCount;
                }

                SupplementaryReadData suppData = SupplementaryReadData.from(read);

                if((suppData.Chromosome.equals(chromosome)) && currentPartition.containsPosition(suppData.Position))
                {
                    expectedRange.setStart(min(expectedRange.start(), suppData.Position));
                    expectedRange.setEnd(max(expectedRange.end(), suppData.Position));
                }
                else
                {
                    chrPartitions.add(formChromosomePartition(suppData.Chromosome, suppData.Position, partitionSize));
                }
            }
        }

        IsComplete = (expectedNonSuppCount == nonSuppCount) && (expectedSuppCount == suppCount);
        ChrPartitions = !chrPartitions.isEmpty() ? chrPartitions : null;
        CurrentRange = currentRange;
        ExpectedRange = expectedRange;
    }

    public static int maxPositionStart(final SAMRecord read)
    {
        int maxPositionStart = read.getAlignmentStart();
        String chromosome = read.getContig();

        if(read.getReadPairedFlag() && !read.getMateUnmappedFlag() && read.getMateReferenceName().equals(chromosome))
            maxPositionStart = max(maxPositionStart, read.getMateAlignmentStart());

        if(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
        {
            SupplementaryReadData suppData = SupplementaryReadData.from(read);

            if((suppData.Chromosome.equals(chromosome)))
                maxPositionStart = max(maxPositionStart, suppData.Position);
        }

        return maxPositionStart;
    }

    public String toString()
    {
        return format("range(%s) expected(%s) complete(%s) dup(%s) remotePartitions(%d) id(%s)",
                CurrentRange, ExpectedRange, IsComplete, Status, ChrPartitions != null ? ChrPartitions.size() : 0, ReadId);
    }
}
