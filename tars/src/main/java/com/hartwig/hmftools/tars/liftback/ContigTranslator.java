package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipped;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightSoftClipped;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.tars.common.ContigEntry;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

// translates a transcript-contig alignment (pos + cigar) to genome coordinates, re-inserting introns as N gaps.
public final class ContigTranslator
{
    public static ContigTranslateResult translate(final ContigEntry contig, final int contigPos, final Cigar contigCigar)
    {
        List<BaseRegion> spans = contig.exonSpans();
        if(spans.isEmpty() || contigPos > contig.contigEnd())
        {
            return null;
        }

        Cigar clampedCigar = clampToContigBounds(contig, contigPos, contigCigar);
        if(clampedCigar == null)
        {
            return null;
        }

        int clampedContigPos = Math.max(contigPos, contig.contigStart());
        int localPos = clampedContigPos - contig.contigStart() + 1;
        SpanLocation start = locateStartSpan(spans, localPos);
        if(start == null)
        {
            return null; // read starts past the end of the contig
        }

        return walkCigarToGenome(contig, spans, start, clampedCigar, hasSoftClipAtExonBoundary(contig, contigPos, contigCigar));
    }

    // bwa-mem2 can anchor into the alt-contig spacer-N region on either side; convert the overhanging M
    // bases to soft-clip so the read still lifts. Returns null when the overhang can't be absorbed.
    private static Cigar clampToContigBounds(final ContigEntry contig, int contigPos, Cigar cigar)
    {
        if(contigPos < contig.contigStart())
        {
            int overhang = contig.contigStart() - contigPos;
            cigar = clampLeadingMToSoftClip(cigar, overhang);
            if(cigar == null)
            {
                return null;
            }
            contigPos = contig.contigStart();
        }

        int readEnd = contigPos + cigar.getReferenceLength() - 1;
        if(readEnd > contig.contigEnd())
        {
            int overhang = readEnd - contig.contigEnd();
            cigar = clampTrailingMToSoftClip(cigar, overhang);
            if(cigar == null)
            {
                return null;
            }
        }

        return cigar;
    }

    // localPos is a 1-based offset into the concatenated exons; returns the containing span and its genomic position.
    private static SpanLocation locateStartSpan(final List<BaseRegion> spans, final int localPos)
    {
        int exonLengthSoFar = 0;
        for(int i = 0; i < spans.size(); ++i)
        {
            BaseRegion span = spans.get(i);
            if(localPos <= exonLengthSoFar + span.baseLength())
            {
                int genomicPos = span.start() + (localPos - exonLengthSoFar - 1);
                return new SpanLocation(i, genomicPos);
            }
            exonLengthSoFar += span.baseLength();
        }
        return null;
    }

    // Walks the clamped contig-space CIGAR, emitting an N at every exon boundary crossed to produce a
    // genome-spaced CIGAR. Returns null when the read extends past the last exon.
    private static ContigTranslateResult walkCigarToGenome(
            final ContigEntry contig, final List<BaseRegion> spans, final SpanLocation start, final Cigar cigar,
            final boolean softClipAtExonBoundary)
    {
        int genomicStart = start.genomicPos();
        int currentSpanIndex = start.spanIndex();
        int currentGenomicPos = start.genomicPos();

        List<CigarElement> outElements = new ArrayList<>();

        for(CigarElement element : cigar.getCigarElements())
        {
            CigarOperator op = element.getOperator();
            int remaining = element.getLength();

            if(!op.consumesReferenceBases())
            {
                outElements.add(new CigarElement(remaining, op));
                continue;
            }

            while(remaining > 0)
            {
                BaseRegion currentSpan = spans.get(currentSpanIndex);

                if(currentGenomicPos > currentSpan.end())
                {
                    if(currentSpanIndex + 1 >= spans.size())
                    {
                        return null; // read extends past last span
                    }

                    BaseRegion nextSpan = spans.get(currentSpanIndex + 1);
                    int intronStart = currentSpan.end() + 1;
                    int intronEnd = nextSpan.start() - 1;

                    if(intronEnd >= intronStart)
                    {
                        outElements.add(new CigarElement(intronEnd - intronStart + 1, CigarOperator.N));
                    }

                    ++currentSpanIndex;
                    currentSpan = nextSpan;
                    currentGenomicPos = currentSpan.start();
                }

                int remainInSpan = currentSpan.end() - currentGenomicPos + 1;
                int take = Math.min(remaining, remainInSpan);

                outElements.add(new CigarElement(take, op));
                currentGenomicPos += take;
                remaining -= take;
            }
        }

        return new ContigTranslateResult(
                contig.chromosome(), genomicStart,
                new Cigar(mergeAdjacentSameOp(dropZeroLength(outElements))),
                softClipAtExonBoundary);
    }

    static List<CigarElement> dropZeroLength(final List<CigarElement> elements)
    {
        List<CigarElement> result = new ArrayList<>(elements.size());

        for(CigarElement element : elements)
        {
            if(element.getLength() > 0)
                result.add(element);
        }

        return result;
    }

    static List<CigarElement> mergeAdjacentSameOp(final List<CigarElement> elements)
    {
        List<CigarElement> merged = new ArrayList<>(elements.size());

        for(CigarElement element : elements)
        {
            if(!merged.isEmpty() && merged.get(merged.size() - 1).getOperator() == element.getOperator())
            {
                CigarElement prev = merged.remove(merged.size() - 1);
                merged.add(new CigarElement(prev.getLength() + element.getLength(), element.getOperator()));
            }
            else
            {
                merged.add(element);
            }
        }

        return merged;
    }

    // exon span index + genomic position where an alignment begins.
    private record SpanLocation(int spanIndex, int genomicPos)
    {
    }

    // true if a leading or trailing S abuts an interior exon boundary (not the outermost end of the first/last exon)
    private static boolean hasSoftClipAtExonBoundary(
            final ContigEntry contig, final int contigPos, final Cigar contigCigar)
    {
        List<BaseRegion> spans = contig.exonSpans();
        if(spans.size() < 2 || contigCigar.isEmpty())
        {
            return false;
        }

        boolean leadingSoftClip = leftSoftClipped(contigCigar);
        boolean trailingSoftClip = rightSoftClipped(contigCigar);
        if(!leadingSoftClip && !trailingSoftClip)
        {
            return false;
        }

        int localPos = contigPos - contig.contigStart() + 1;
        int endLocalPos = localPos + contigCigar.getReferenceLength() - 1;

        int exonLengthSoFar = 0;
        for(int i = 0; i < spans.size() - 1; ++i)
        {
            exonLengthSoFar += spans.get(i).baseLength();
            if(leadingSoftClip && localPos == exonLengthSoFar + 1)
            {
                return true;
            }
            if(trailingSoftClip && endLocalPos == exonLengthSoFar)
            {
                return true;
            }
        }
        return false;
    }

    // converts the first `overhang` M bases to soft-clip, merging any existing leading S; null if not possible
    private static Cigar clampLeadingMToSoftClip(final Cigar cigar, final int overhang)
    {
        List<CigarElement> elements = cigar.getCigarElements();

        int existingLeadingSoftClip = leftSoftClipLength(cigar);
        int alignedIndex = existingLeadingSoftClip > 0 ? 1 : 0;

        CigarElement aligned = alignedIndex < elements.size() ? elements.get(alignedIndex) : null;
        if(aligned == null || aligned.getOperator() != CigarOperator.M || aligned.getLength() <= overhang)
        {
            return null;
        }

        List<CigarElement> out = new ArrayList<>(elements.size());
        out.add(new CigarElement(existingLeadingSoftClip + overhang, CigarOperator.S));
        out.add(new CigarElement(aligned.getLength() - overhang, CigarOperator.M));
        out.addAll(elements.subList(alignedIndex + 1, elements.size()));
        return new Cigar(out);
    }

    // trailing-edge mirror of clampLeadingMToSoftClip: reverse the elements, clamp the (now-leading) edge, reverse back.
    private static Cigar clampTrailingMToSoftClip(final Cigar cigar, final int overhang)
    {
        List<CigarElement> reversed = new ArrayList<>(cigar.getCigarElements());
        Collections.reverse(reversed);

        Cigar clamped = clampLeadingMToSoftClip(new Cigar(reversed), overhang);
        if(clamped == null)
        {
            return null;
        }

        List<CigarElement> out = new ArrayList<>(clamped.getCigarElements());
        Collections.reverse(out);
        return new Cigar(out);
    }
}
