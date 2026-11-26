package com.hartwig.hmftools.tars.liftback.emit;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.tars.liftback.translate.ContigTranslator;

import htsjdk.samtools.SAMRecord;

// Rebuilds reciprocal SA tags from the final records in one name group. The caller applies placements and refreshes NM
// first, so an entry cannot drift from the record that is emitted.
public final class SaTagRewriter
{
    // Admission check only: a supplementary with no liftable SA partner is orphaned. The emitted SA is rebuilt from the
    // final name group, not from these pre-lift coordinates.
    static boolean hasLiftableAlignment(final String saTag, final ContigTranslator translator)
    {
        List<SupplementaryReadData> alignments;
        try
        {
            alignments = SupplementaryReadData.extractAlignments(saTag);
        }
        catch(RuntimeException e)
        {
            alignments = null;
        }

        if(alignments == null)
        {
            throw new IllegalStateException(String.format("malformed SA(%s)", saTag));
        }

        for(SupplementaryReadData alignment : alignments)
        {
            if(translator.liftAlignment(
                    alignment.Chromosome, alignment.Position, alignment.Cigar,
                    alignment.NM, alignment.isForwardOrient()) != null)
            {
                return true;
            }
        }
        return false;
    }

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
