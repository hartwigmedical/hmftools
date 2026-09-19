package com.hartwig.hmftools.virusdetect;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightClipLength;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.virusdetect.VirusConstants.ORIGIN_CLIP_TOLERANCE;

import java.util.Comparator;
import java.util.List;

import htsjdk.samtools.SAMRecord;

// One viral alignment of a read with just the data required for our analysis.
public record ViralAlignment(
        String readName,
        ViralContig contig,
        int alignmentStart,
        int alignmentEnd,
        int leftClip,
        int rightClip,
        int alignerScore,
        // NM + clipped bases: how poorly the contig explains the read
        int divergence,
        List<AlignedInterval> alignedIntervals)
{
    public ViralAlignment
    {
        if(readName.isEmpty())
        {
            throw new IllegalArgumentException("read name is empty");
        }
        if(alignmentStart < 1)
        {
            throw new IllegalArgumentException("invalid alignment start: " + alignmentStart);
        }
        if(alignmentEnd < alignmentStart)
        {
            throw new IllegalArgumentException("invalid alignment span: " + alignmentStart + "-" + alignmentEnd);
        }
        if(leftClip < 0 || rightClip < 0)
        {
            throw new IllegalArgumentException("invalid clips: left " + leftClip + ", right " + rightClip);
        }
        if(alignerScore < 0)
        {
            throw new IllegalArgumentException("invalid aligner score: " + alignerScore);
        }
        if(divergence < leftClip + rightClip)
        {
            throw new IllegalArgumentException("divergence below clipped bases: " + divergence);
        }
        if(alignedIntervals.isEmpty())
        {
            throw new IllegalArgumentException("no aligned intervals");
        }
    }

    public static ViralAlignment from(SAMRecord record, ViralReference reference)
    {
        int leftClip = leftClipLength(record.getCigar());
        int rightClip = rightClipLength(record.getCigar());

        List<AlignedInterval> intervals = record.getAlignmentBlocks().stream()
                .map(block -> new AlignedInterval(block.getReferenceStart(), block.getLength()))
                .toList();

        int editDistance = requiredTag(record, NUM_MUTATONS_ATTRIBUTE, "edit distance (NM)");
        int alignerScore = requiredTag(record, ALIGNMENT_SCORE_ATTRIBUTE, "alignment score");

        return new ViralAlignment(
                record.getReadName(), reference.contig(record.getReferenceName()), record.getAlignmentStart(),
                record.getAlignmentEnd(), leftClip, rightClip, alignerScore, editDistance + leftClip + rightClip, intervals);
    }

    // The clipped bases project past a contig end, so the read straddles the circular genome's linearization origin:
    // those bases wrap to the other end and cannot align linearly. An artifact, not real divergence.
    public boolean clipsOverContigEnd()
    {
        int contigLength = contig.length();
        return alignmentStart - leftClip < 1 - ORIGIN_CLIP_TOLERANCE || alignmentEnd + rightClip > contigLength + ORIGIN_CLIP_TOLERANCE;
    }

    private static int requiredTag(SAMRecord record, String tag, String description)
    {
        Integer value = record.getIntegerAttribute(tag);
        if(value == null)
        {
            throw new IllegalStateException(
                    String.format("aligned read missing %s on contig: %s", description, record.getReferenceName()));
        }
        return value;
    }

    // Sort primarily by divergence because we are mostly interested in per-base difference.
    // Aligner score is not the best fit.
    public static final Comparator<ViralAlignment> BEST_FIT_FIRST = Comparator
            .comparingInt(ViralAlignment::divergence)
            .thenComparing(ViralAlignment::alignerScore, Comparator.reverseOrder())
            .thenComparingInt(ViralAlignment::alignmentStart);

    // TODO: should be separate file
    // A contiguous run of reference bases covered by the alignment (a CIGAR M/=/X block); 1-based reference start.
    public record AlignedInterval(
            int referenceStart,
            int length
    )
    {
        public AlignedInterval
        {
            if(referenceStart < 1)
            {
                throw new IllegalArgumentException("invalid reference start: " + referenceStart);
            }
            if(length < 1)
            {
                throw new IllegalArgumentException("invalid interval length: " + length);
            }
        }
    }
}
