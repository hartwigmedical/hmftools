package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_NEG_STRAND;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;

import java.util.List;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import htsjdk.samtools.SAMRecord;

// Rebuilds reciprocal SA tags from the final records in one name group. The caller applies placements and refreshes NM
// first, so an entry cannot drift from the record that is emitted.
public final class SaTagRewriter
{
    private SaTagRewriter()
    {
    }

    // Admission check only: a supplementary with no liftable SA partner is orphaned. The emitted SA is rebuilt from the
    // final name group, not from these pre-lift coordinates.
    static boolean hasLiftableAlignment(final String saTag, final ContigTranslator translator)
    {
        List<SupplementaryReadData> alignments = SupplementaryReadData.extractAlignments(saTag);
        if(alignments == null)
        {
            return false;
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
        StringBuilder saTag = new StringBuilder();

        for(int i = 0; i < records.size(); ++i)
        {
            if(i == recordIndex || !willEmit[i])
            {
                continue;
            }

            SAMRecord alignment = records.get(i);
            if(alignment.getReadUnmappedFlag())
            {
                continue;
            }

            boolean include = targetIsSupplementary
                    ? alignment == primary || alignment.getSupplementaryAlignmentFlag()
                    : alignment.getSupplementaryAlignmentFlag();
            if(include)
            {
                saTag.append(asSaEntry(alignment));
            }
        }

        return saTag.length() > 0 ? saTag.toString() : null;
    }

    private static String asSaEntry(final SAMRecord record)
    {
        Integer nm = record.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE);
        return record.getReferenceName() + ',' + record.getAlignmentStart()
                + ',' + (record.getReadNegativeStrandFlag() ? SUPP_NEG_STRAND : SUPP_POS_STRAND)
                + ',' + record.getCigarString().replace('H', 'S') + ',' + record.getMappingQuality()
                + ',' + (nm != null ? nm : 0) + ';';
    }
}
