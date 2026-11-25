package com.hartwig.hmftools.tars.liftback.translate;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.XA_ATTRIBUTE;
import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.tars.common.ContigEntry;
import com.hartwig.hmftools.tars.liftback.AlignmentKey;
import com.hartwig.hmftools.tars.liftback.LiftedAlignment;
import com.hartwig.hmftools.tars.liftback.LiftedRecord;
import com.hartwig.hmftools.tars.liftback.overhang.OverhangGate;

import htsjdk.samtools.SAMRecord;

public final class AlignmentLifter
{
    private final ContigTranslator mTranslator;

    public AlignmentLifter(final List<ContigEntry> entries)
    {
        mTranslator = new ContigTranslator(entries);
    }

    public ContigTranslator translator()
    {
        return mTranslator;
    }

    public LiftedRecord liftPrimary(final SAMRecord record, final OverhangGate overhangGate)
    {
        if(record.getReadUnmappedFlag())
        {
            return LiftedRecord.unmapped("");
        }

        LiftedAlignment self = liftSelf(record);
        if(self == null)
        {
            return LiftedRecord.unmapped("primary_translate_failed");
        }

        List<LiftedAlignment> alternatives =
                mTranslator.liftXaAlignments(record.getStringAttribute(XA_ATTRIBUTE));
        List<LiftedAlignment> alignments = new ArrayList<>(1 + alternatives.size());
        alignments.add(self);
        alignments.addAll(alternatives);
        overhangGate.gatePlacements(alignments, record);

        return new LiftedRecord(record.getMappingQuality(), 0, "", 0, alignments);
    }

    public LiftedRecord liftSupplementary(final SAMRecord record, final OverhangGate overhangGate)
    {
        LiftedAlignment self = liftSelf(record);
        if(self == null)
        {
            return LiftedRecord.unmapped("supp_translate_failed");
        }

        List<LiftedAlignment> alignments = new ArrayList<>();
        Set<AlignmentKey> seen = new HashSet<>();
        alignments.add(self);
        seen.add(self.key());
        for(LiftedAlignment alignment : mTranslator.liftXaAlignments(record.getStringAttribute(XA_ATTRIBUTE)))
        {
            if(seen.add(alignment.key()))
            {
                alignments.add(alignment);
            }
        }
        overhangGate.gatePlacements(alignments, record);

        return new LiftedRecord(record.getMappingQuality(), 1, "", 0, alignments);
    }

    private LiftedAlignment liftSelf(final SAMRecord record)
    {
        Integer mismatches = record.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE);
        LiftedAlignment lifted = mTranslator.liftAlignment(
                record.getReferenceName(), record.getAlignmentStart(), record.getCigarString(),
                mismatches != null ? mismatches : 0, !record.getReadNegativeStrandFlag());
        if(lifted == null)
        {
            logLiftFailure(record);
        }
        return lifted;
    }

    private void logLiftFailure(final SAMRecord record)
    {
        String contig = record.getReferenceName();
        int pos = record.getAlignmentStart();
        int readEnd = pos + record.getCigar().getReferenceLength() - 1;
        String role = record.getSupplementaryAlignmentFlag() ? "supp" : "primary";
        ContigEntry segment = mTranslator.findSegment(contig, pos);

        if(segment == null)
        {
            TARS_LOGGER.debug(
                    "lift failed {}: {}:{}-{} {} - contig has no segments",
                    role, contig, pos, readEnd, record.getCigarString());
            return;
        }

        // findSegment returns the nearest segment for spacer positions.
        if(pos < segment.contigStart() || pos > segment.contigEnd())
        {
            return;
        }

        TARS_LOGGER.debug(
                "lift failed {}: {}:{}-{} {} - segment {} [{}-{}] exons({})",
                role, contig, pos, readEnd, record.getCigarString(),
                segment.transName(), segment.contigStart(), segment.contigEnd(), segment.exonSpans().size());
    }
}
