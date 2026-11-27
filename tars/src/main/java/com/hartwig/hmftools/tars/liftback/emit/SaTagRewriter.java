package com.hartwig.hmftools.tars.liftback.emit;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import htsjdk.samtools.SAMRecord;

// Rebuilds reciprocal SA tags from the final records in one name group. The caller applies placements and refreshes NM
// first, so an entry cannot drift from the record that is emitted.
public final class SaTagRewriter
{
    static String buildForRecord(
            final List<SAMRecord> records, final boolean[] willEmit, final SAMRecord primary, final int recordIndex)
    {
        SAMRecord target = records.get(recordIndex);
        boolean targetIsSupplementary = target.getSupplementaryAlignmentFlag();
        List<SupplementaryReadData> alignments = new ArrayList<>();

        for(int i = 0; i < records.size(); ++i)
        {
            if(i == recordIndex || !willEmit[i])
            {
                continue;
            }

            SAMRecord alignment = records.get(i);
            boolean include = targetIsSupplementary
                    ? alignment == primary || alignment.getSupplementaryAlignmentFlag()
                    : alignment.getSupplementaryAlignmentFlag();
            if(include)
            {
                alignments.add(SupplementaryReadData.createSupplementaryReadDataForRecord(alignment));
            }
        }

        return alignments.isEmpty() ? null : SupplementaryReadData.alignmentsToSamTag(alignments);
    }
}
