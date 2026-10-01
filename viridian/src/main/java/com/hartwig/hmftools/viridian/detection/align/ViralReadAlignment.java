package com.hartwig.hmftools.viridian.detection.align;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightClipLength;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRAL_CONTIG_ORIGIN_CLIP_TOLERANCE;

import java.util.Comparator;
import java.util.List;

import com.hartwig.hmftools.viridian.detection.common.ReadId;
import com.hartwig.hmftools.viridian.reference.ViralContig;
import com.hartwig.hmftools.viridian.reference.VirusReference;

import htsjdk.samtools.SAMRecord;

// One alignment of a read to a virus genome, with just the data required for our analysis, in a nice format.
public record ViralReadAlignment(
        ReadId readId,
        ViralContig contig,
        int alignmentStart,
        int alignmentEnd,
        String cigar,
        int leftClip,
        int rightClip,
        int alignerScore,
        // NM + clipped bases: how poorly the contig explains the read
        int divergence,
        List<AlignedInterval> alignedIntervals)
{
    public ViralReadAlignment
    {
        if(alignmentStart < 1)
        {
            throw new IllegalArgumentException("Invalid alignmentStart: " + alignmentStart);
        }
        if(alignmentEnd < alignmentStart)
        {
            throw new IllegalArgumentException("Invalid alignment span: " + alignmentStart + "-" + alignmentEnd);
        }
        if(leftClip < 0 || rightClip < 0)
        {
            throw new IllegalArgumentException("Invalid clips: left " + leftClip + ", right " + rightClip);
        }
        if(alignerScore < 0)
        {
            throw new IllegalArgumentException("Invalid alignerScore: " + alignerScore);
        }
        if(divergence < leftClip + rightClip)
        {
            throw new IllegalArgumentException("Invalid divergence: " + divergence);
        }
        if(cigar.isEmpty())
        {
            throw new IllegalArgumentException("Invalid cigar");
        }
        if(alignedIntervals.isEmpty())
        {
            throw new IllegalArgumentException("Invalid alignedIntervals");
        }
    }

    public static ViralReadAlignment from(SAMRecord record, VirusReference reference)
    {
        int leftClip = leftClipLength(record.getCigar());
        int rightClip = rightClipLength(record.getCigar());

        List<AlignedInterval> intervals = record.getAlignmentBlocks().stream()
                .map(block -> new AlignedInterval(block.getReferenceStart(), block.getLength()))
                .toList();

        int editDistance = requiredTag(record, NUM_MUTATONS_ATTRIBUTE, "edit distance (NM)");
        int alignerScore = requiredTag(record, ALIGNMENT_SCORE_ATTRIBUTE, "alignment score");

        return new ViralReadAlignment(
                ReadId.parse(record.getReadName()), reference.contig(record.getReferenceName()), record.getAlignmentStart(),
                record.getAlignmentEnd(), record.getCigarString(), leftClip, rightClip, alignerScore, editDistance + leftClip + rightClip,
                intervals);
    }

    // The clipped bases project past a contig end, so the read straddles the circular genome's linearization origin:
    // those bases wrap to the other end and cannot align linearly. An artifact, not real divergence.
    public boolean clipsOverContigEnd()
    {
        int contigLength = contig.length();
        return alignmentStart - leftClip < 1 - VIRAL_CONTIG_ORIGIN_CLIP_TOLERANCE
                || alignmentEnd + rightClip > contigLength + VIRAL_CONTIG_ORIGIN_CLIP_TOLERANCE;
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
    public static final Comparator<ViralReadAlignment> BEST_FIT_FIRST = Comparator
            .comparingInt(ViralReadAlignment::divergence)
            .thenComparing(ViralReadAlignment::alignerScore, Comparator.reverseOrder())
            .thenComparingInt(ViralReadAlignment::alignmentStart);
}
