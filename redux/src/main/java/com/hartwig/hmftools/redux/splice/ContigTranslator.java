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

        // bwa-mem2 can anchor into the alt-contig spacer-N region on either side; reclaim overhanging M
        // bases as soft-clip so the read still lifts. Returns null when the overhang can't be absorbed.
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

        // walk contig-space CIGAR, emitting N at each exon boundary to produce a genome-spaced CIGAR
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

    // A D straddling an exon boundary lifts as xD nN yD — artifact of tx-FASTA off-by-N at the junction.
    // Absorb small flanking Ds into the N (preserves read and ref span); keep larger ones as real deletions.
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

            if(!result.isEmpty() && isAbsorbableDeletion(result.get(result.size() - 1)))
                splicedLength += result.remove(result.size() - 1).getLength();

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

    // Drops fabricated sub-threshold junction anchors at both ends. A tx-contig walk that dribbles past
    // an exon boundary leaves a tiny "yM nN..." or "...nN yM zS" anchor; STAR rejects junctions whose
    // flanking block is shorter than alignSJDBoverhangMin. Two floors: bareAnchorFloor for a read
    // genuinely starting inside an exon (no adjacent S, kept down to 1bp); softclipAnchorFloor for
    // an anchor sitting next to a softclip (bwa over-ran the boundary, so the implied junction is
    // unsupported — roll the yM into the softclip). Leading trim advances start past dropped yM+intron.
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
        // without a preceding M the trimmed cigar would start with S/N — refuse
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

    // "[zS]? yM nN <...>M<...>" — leading S optional (absent when source starts with M, no clip to merge)
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

    // true if a leading or trailing S abuts an interior exon boundary (not the outermost end of the first/last exon)
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

        final int localPos = contigPos - contig.altStart() + 1;
        final int endLocalPos = localPos + contigCigar.getReferenceLength() - 1;

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

    // converts the first `overhang` M bases to soft-clip, merging any existing leading S; null if not possible
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

    // trailing-edge mirror of clampLeadingMToSoftClip
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
