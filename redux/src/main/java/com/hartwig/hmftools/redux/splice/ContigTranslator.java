package com.hartwig.hmftools.redux.splice;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public final class ContigTranslator
{
    // tx-contig D ops up to this length that end up adjacent to a lift-emitted N are absorbed into the N
    // (see collapseSpliceFlankingDeletions). Larger D values are kept as real deletions.
    static final int SPLICE_FLANKING_DELETION_MAX_BP = 5;

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

        return new TranslationResult(
                contig.chromosome(), genomicStart,
                new Cigar(collapseSpliceFlankingDeletions(outElements)),
                impliedIntrons);
    }

    // A pre-lift D op that straddles an exon boundary gets emitted by translate() as xD nN yD (lead D
    // up to the exon end, splice N, trail D into the next exon). This is the signature of a tx-FASTA
    // off-by-N at the junction; the D consumes no read bases and the N consumes no read bases either,
    // so absorbing xD..yD into the N preserves read and reference span. Only absorbed when each D's
    // length is <= SPLICE_FLANKING_DELETION_MAX_BP — larger Ds are kept as real deletions.
    static List<CigarElement> collapseSpliceFlankingDeletions(final List<CigarElement> elements)
    {
        if(elements.size() < 3)
            return elements;

        List<CigarElement> result = new ArrayList<>(elements.size());
        for(int i = 0; i < elements.size(); ++i)
        {
            final CigarElement element = elements.get(i);

            if(element.getOperator() != CigarOperator.N)
            {
                result.add(element);
                continue;
            }

            int splicedLength = element.getLength();

            // absorb a leading absorbable D already emitted into result
            if(!result.isEmpty() && isAbsorbableDeletion(result.get(result.size() - 1)))
                splicedLength += result.remove(result.size() - 1).getLength();

            // absorb a trailing absorbable D that immediately follows in the input
            while(i + 1 < elements.size() && isAbsorbableDeletion(elements.get(i + 1)))
                splicedLength += elements.get(++i).getLength();

            result.add(new CigarElement(splicedLength, CigarOperator.N));
        }
        return result;
    }

    private static boolean isAbsorbableDeletion(final CigarElement element)
    {
        return element.getOperator() == CigarOperator.D && element.getLength() <= SPLICE_FLANKING_DELETION_MAX_BP;
    }

    public static final class MicroAnchorResult
    {
        public final Cigar AdjustedCigar;
        public final int StartShift;     // bases the alignment start advanced (from a leading trim)

        MicroAnchorResult(final Cigar cigar, final int startShift)
        {
            AdjustedCigar = cigar;
            StartShift = startShift;
        }
    }

    // Drops fabricated sub-threshold junction anchors at both read ends. A tx-contig walk that
    // dribbles a few bases past an exon boundary leaves a tiny terminal "yM nN ..." (leading) or
    // "... nN yM zS" (trailing) anchor. STAR rejects an annotated junction whose flanking exon
    // block is shorter than alignSJDBoverhangMin (stitchWindowAligns.cpp), so we do the same:
    // when y < minAnchorBp the yM is rolled into the adjacent softclip and the junction dropped.
    // A leading trim advances the alignment start past the dropped yM + intron (reported via
    // StartShift). A "<...>N yM" tail with no softclip is left intact — there the short exon is
    // the read's genuine end.
    // Two floors. bareAnchorFloor applies to a tiny anchor at the read's true terminus (a bare leading
    // yM nN... with no leading softclip): that yM is where the read genuinely starts inside an exon, so
    // it is kept down to 1bp. softclipAnchorFloor applies to a tiny anchor sitting NEXT TO a softclip
    // (...nN yM zS trailing, or zS yM nN... leading): there the read did not actually span the junction -
    // bwa over-ran the exon boundary by a few bases and clipped the rest - so the implied junction is
    // unsupported and the tiny anchor is rolled into the softclip below this floor.
    public static MicroAnchorResult trimMicroAnchors(
            final Cigar cigar, final int bareAnchorFloor, final int softclipAnchorFloor)
    {
        List<CigarElement> elements = new ArrayList<>(cigar.getCigarElements());
        elements = trimTrailingAnchor(elements, softclipAnchorFloor);   // trailing trim is always softclip-adjacent
        int startShift = 0;

        final boolean hasLeadingS = !elements.isEmpty() && elements.get(0).getOperator() == CigarOperator.S;
        final int leadingFloor = hasLeadingS ? softclipAnchorFloor : bareAnchorFloor;
        if(leadingAnchorTrimmable(elements, leadingFloor))
        {
            final int anchorIdx = hasLeadingS ? 1 : 0;
            final int existingS = hasLeadingS ? elements.get(0).getLength() : 0;
            final int tinyAnchor = elements.get(anchorIdx).getLength();
            final int intron = elements.get(anchorIdx + 1).getLength();
            startShift = tinyAnchor + intron;
            final List<CigarElement> trimmed = new ArrayList<>(elements.size() - anchorIdx - 1);
            trimmed.add(new CigarElement(existingS + tinyAnchor, CigarOperator.S));
            for(int i = anchorIdx + 2; i < elements.size(); ++i)
                trimmed.add(elements.get(i));
            elements = trimmed;
        }

        return new MicroAnchorResult(new Cigar(elements), startShift);
    }

    // "<...>M nN yM zS" -> "<...>M (y+z)S" when y < minAnchorBp.
    public static Cigar trimTrailingMicroAnchor(final Cigar cigar, final int minAnchorBp)
    {
        return new Cigar(trimTrailingAnchor(new ArrayList<>(cigar.getCigarElements()), minAnchorBp));
    }

    private static List<CigarElement> trimTrailingAnchor(final List<CigarElement> elements, final int minAnchorBp)
    {
        final int last = elements.size() - 1;
        if(last < 2)
            return elements;
        if(elements.get(last).getOperator() != CigarOperator.S)
            return elements;
        final CigarElement tailMatch = elements.get(last - 1);
        if(tailMatch.getOperator() != CigarOperator.M || tailMatch.getLength() >= minAnchorBp)
            return elements;
        if(elements.get(last - 2).getOperator() != CigarOperator.N)
            return elements;
        // refuse if there's no preceding M anchor — would leave the cigar starting with S/N.
        boolean hasAnchorBeforeIntron = false;
        for(int i = 0; i < last - 2; ++i)
        {
            if(elements.get(i).getOperator() == CigarOperator.M)
            {
                hasAnchorBeforeIntron = true;
                break;
            }
        }
        if(!hasAnchorBeforeIntron)
            return elements;

        final List<CigarElement> result = new ArrayList<>(last - 1);
        for(int i = 0; i < last - 2; ++i)
            result.add(elements.get(i));
        result.add(new CigarElement(tailMatch.getLength() + elements.get(last).getLength(), CigarOperator.S));
        return result;
    }

    // matches "[zS]? yM nN <...>M<...>" — optional leading softclip, tiny matched anchor under
    // threshold, intron, then a real matched anchor that survives the trim. The leading S is
    // optional because a tx-contig source cigar that begins with M (no clip) lifts to a cigar
    // starting directly with the tiny anchor — same fabrication pattern, no leading S to merge into.
    private static boolean leadingAnchorTrimmable(final List<CigarElement> elements, final int minAnchorBp)
    {
        if(elements.isEmpty())
            return false;
        final boolean hasLeadingS = elements.get(0).getOperator() == CigarOperator.S;
        final int anchorIdx = hasLeadingS ? 1 : 0;
        if(anchorIdx + 2 >= elements.size())
            return false;
        final CigarElement tinyAnchor = elements.get(anchorIdx);
        if(tinyAnchor.getOperator() != CigarOperator.M || tinyAnchor.getLength() >= minAnchorBp)
            return false;
        if(elements.get(anchorIdx + 1).getOperator() != CigarOperator.N)
            return false;
        for(int i = anchorIdx + 2; i < elements.size(); ++i)
        {
            if(elements.get(i).getOperator() == CigarOperator.M)
                return true;
        }
        return false;
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
