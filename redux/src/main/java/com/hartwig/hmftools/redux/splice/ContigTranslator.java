package com.hartwig.hmftools.redux.splice;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public final class ContigTranslator
{
    public static TranslationResult translate(final ContigEntry contig, int contigPos, final Cigar contigCigar)
    {
        final List<BaseRegion> spans = contig.exonSpans();
        if(spans.isEmpty() || contigPos > contig.altEnd())
            return null;

        // bwa-mem2 sometimes anchors a read 1-N bases before this transcript's altStart (i.e. into the
        // upstream alt-packing spacer-N region) or extends a few bases past altEnd into the downstream
        // spacer. Reclaim those overhanging M bases as soft-clip so the read still lifts. Falls through
        // to null when the overhang lands on a non-M op or the M is too short to absorb the shift.
        Cigar workingCigar = contigCigar;
        if(contigPos < contig.altStart())
        {
            int overhang = contig.altStart() - contigPos;
            workingCigar = clampLeadingMToSoftClip(workingCigar, overhang);
            if(workingCigar == null)
                return null;
            contigPos = contig.altStart();
        }

        int readEnd = contigPos + workingCigar.getReferenceLength() - 1;
        if(readEnd > contig.altEnd())
        {
            int overhang = readEnd - contig.altEnd();
            workingCigar = clampTrailingMToSoftClip(workingCigar, overhang);
            if(workingCigar == null)
                return null;
        }

        int localPos = contigPos - contig.altStart() + 1;

        // locate the span that contains the read's leftmost mapped contig position
        int currentSpanIndex = -1;
        int currentGenomicPos = -1;
        int exonLengthSoFar = 0;
        for(int i = 0; i < spans.size(); ++i)
        {
            BaseRegion span = spans.get(i);
            if(localPos <= exonLengthSoFar + span.baseLength())
            {
                currentSpanIndex = i;
                currentGenomicPos = span.start() + (localPos - exonLengthSoFar - 1);
                break;
            }
            exonLengthSoFar += span.baseLength();
        }

        if(currentSpanIndex < 0)
            return null; // read starts past the end of the contig

        int genomicStart = currentGenomicPos;
        List<CigarElement> outElements = new ArrayList<>();
        List<BaseRegion> impliedIntrons = new ArrayList<>();

        // walk the contig-space CIGAR left-to-right; for each reference-consuming chunk, step through the
        // exon spans, copying the operator out and emitting an N whenever we cross an exon boundary so the
        // rebuilt CIGAR is genome-spaced. Non-reference ops (I, S) pass through unchanged.
        for(CigarElement element : workingCigar.getCigarElements())
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

                // if we've stepped past the current span, insert an intron and advance into the next span
                if(currentGenomicPos > currentSpan.end())
                {
                    if(currentSpanIndex + 1 >= spans.size())
                        return null; // read extends past last span

                    BaseRegion nextSpan = spans.get(currentSpanIndex + 1);
                    int intronStart = currentSpan.end() + 1;
                    int intronEnd = nextSpan.start() - 1;

                    outElements.add(new CigarElement(intronEnd - intronStart + 1, CigarOperator.N));
                    impliedIntrons.add(new BaseRegion(intronStart, intronEnd));

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

        return new TranslationResult(contig.chromosome(), genomicStart, new Cigar(outElements), impliedIntrons);
    }

    // true if the contig CIGAR has a leading or trailing S whose adjacent edge sits at an interior exon
    // boundary of the contig (i.e. a real splice-junction edge, not the outermost end of the first/last exon).
    public static boolean hasSoftClipAtExonBoundary(
            final ContigEntry contig, final int contigPos, final Cigar contigCigar)
    {
        final List<BaseRegion> spans = contig.exonSpans();
        if(spans.size() < 2 || contigCigar.isEmpty())
            return false;

        final List<CigarElement> elements = contigCigar.getCigarElements();
        final boolean leadingSoftClip = elements.get(0).getOperator() == CigarOperator.S;
        final boolean trailingSoftClip = elements.get(elements.size() - 1).getOperator() == CigarOperator.S;
        if(!leadingSoftClip && !trailingSoftClip)
            return false;

        // shift alt-contig position into transcript-local coordinates (1-based)
        final int localPos = contigPos - contig.altStart() + 1;
        final int endLocalPos = localPos + contigCigar.getReferenceLength() - 1;

        // interior exon-end boundary in contig-local coords = cumulative length of exons[0..i] for i < last
        int exonLengthSoFar = 0;
        for(int i = 0; i < spans.size() - 1; ++i)
        {
            exonLengthSoFar += spans.get(i).baseLength();
            if(leadingSoftClip && localPos == exonLengthSoFar + 1)
                return true;
            if(trailingSoftClip && endLocalPos == exonLengthSoFar)
                return true;
        }
        return false;
    }

    // converts the first `overhang` ref-consuming bases of the leading M into soft-clip. Any pre-existing
    // leading S is merged into one larger S. Returns null when the first ref-consuming op is not a long
    // enough M (e.g. starts with an I/D/N, or M is shorter than the overhang).
    private static Cigar clampLeadingMToSoftClip(final Cigar cigar, final int overhang)
    {
        List<CigarElement> elements = cigar.getCigarElements();
        int existingLeadingSoftClip = 0;
        int idx = 0;
        while(idx < elements.size() && elements.get(idx).getOperator() == CigarOperator.S)
        {
            existingLeadingSoftClip += elements.get(idx).getLength();
            ++idx;
        }
        if(idx >= elements.size())
            return null;

        CigarElement first = elements.get(idx);
        if(first.getOperator() != CigarOperator.M || first.getLength() <= overhang)
            return null;

        List<CigarElement> out = new ArrayList<>(elements.size());
        out.add(new CigarElement(existingLeadingSoftClip + overhang, CigarOperator.S));
        out.add(new CigarElement(first.getLength() - overhang, CigarOperator.M));
        for(int i = idx + 1; i < elements.size(); ++i)
            out.add(elements.get(i));
        return new Cigar(out);
    }

    // mirror of clampLeadingMToSoftClip for the trailing edge.
    private static Cigar clampTrailingMToSoftClip(final Cigar cigar, final int overhang)
    {
        List<CigarElement> elements = cigar.getCigarElements();
        int existingTrailingSoftClip = 0;
        int idx = elements.size() - 1;
        while(idx >= 0 && elements.get(idx).getOperator() == CigarOperator.S)
        {
            existingTrailingSoftClip += elements.get(idx).getLength();
            --idx;
        }
        if(idx < 0)
            return null;

        CigarElement last = elements.get(idx);
        if(last.getOperator() != CigarOperator.M || last.getLength() <= overhang)
            return null;

        List<CigarElement> out = new ArrayList<>(elements.size());
        for(int i = 0; i < idx; ++i)
            out.add(elements.get(i));
        out.add(new CigarElement(last.getLength() - overhang, CigarOperator.M));
        out.add(new CigarElement(existingTrailingSoftClip + overhang, CigarOperator.S));
        return new Cigar(out);
    }

    public record TranslationResult(
            String chromosome,
            int genomicStart,
            Cigar genomicCigar,
            List<BaseRegion> impliedIntrons) {}
}
