package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.fusion.FusionRead.convertReads;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.findSplitReadJunction;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.isofox.common.ReadRecord;

public final class FusionTestUtils
{
    public static FusionFragment fromReads(final List<ReadRecord> reads)
    {
        setReadJunctions(reads);
        return new FusionFragment(new FusionReadGroup(reads.get(0).Id, convertReads(reads)));
    }

    public static void setReadJunctions(final List<ReadRecord> reads)
    {
        for(ReadRecord read : reads)
        {
            if(read.containsSplit() && read.junctionPositions() == null)
            {
                // mirror what is done in CRT
                int[] junctionPositions = findSplitReadJunction(read);
                if(junctionPositions != null)
                {
                    read.setJunctionPosition(SE_START, junctionPositions[SE_START]);
                    read.setJunctionPosition(SE_END, junctionPositions[SE_END]);
                }
            }
        }
    }

    public static FusionReadGroup createGroup(final ReadRecord read)
    {
        final List<ReadRecord> reads = Lists.newArrayList(read);
        setReadJunctions(reads);
        List<FusionRead> fusionReads = FusionRead.convertReads(reads);
        return new FusionReadGroup(read.Id, fusionReads);
    }

    public static FusionReadGroup createGroup(final ReadRecord read1, final ReadRecord read2)
    {
        final List<ReadRecord> reads = Lists.newArrayList(read1, read2);
        setReadJunctions(reads);
        List<FusionRead> fusionReads = FusionRead.convertReads(reads);
        return new FusionReadGroup(read1.Id, fusionReads);
    }


}
