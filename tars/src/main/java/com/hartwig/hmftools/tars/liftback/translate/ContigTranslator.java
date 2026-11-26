package com.hartwig.hmftools.tars.liftback.translate;

import static com.hartwig.hmftools.tars.common.TarsCigarUtils.clampLeadingReferenceToSoftClip;
import static com.hartwig.hmftools.tars.common.TarsCigarUtils.clampTrailingReferenceToSoftClip;
import static com.hartwig.hmftools.tars.common.TarsCigarUtils.isMatchOrEqualOp;
import static com.hartwig.hmftools.tars.common.TarsCigarUtils.normalize;
import static com.hartwig.hmftools.tars.common.TarsConstants.ALT_CONTIG_SUFFIX;
import static com.hartwig.hmftools.tars.common.TarsConstants.MAX_MERGED_DELETION_BP;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.tars.common.ContigEntry;
import com.hartwig.hmftools.tars.liftback.AlignmentKey;
import com.hartwig.hmftools.tars.liftback.LiftedAlignment;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

// Translates a transcript-contig alignment (pos + cigar) to genome coordinates, re-inserting introns as N gaps.
// Holds the sidecar segments per alt contig so a position resolves to its owning transcript first; translate() is the
// coordinate walk on its own and needs no index.
public final class ContigTranslator
{
    // sorted by contigStart for binary search back to the owning transcript.
    private final Map<String, List<ContigEntry>> mSegmentsByAltContig;

    public ContigTranslator(final List<ContigEntry> entries)
    {
        mSegmentsByAltContig = new HashMap<>();
        for(ContigEntry entry : entries)
        {
            mSegmentsByAltContig.computeIfAbsent(entry.contigRegion().Chromosome, key -> new ArrayList<>()).add(entry);
        }
        for(List<ContigEntry> segments : mSegmentsByAltContig.values())
        {
            segments.sort(Comparator.comparingInt(entry -> entry.contigRegion().start()));
        }
    }

    public Set<String> contigNames()
    {
        return mSegmentsByAltContig.keySet();
    }

    // Lifts one placement to genomic coordinates. Returns null for an unknown alt contig or a position outside any
    // transcript; ref-genome placements pass through unchanged.
    public LiftedAlignment liftAlignment(
            final String contig, final int pos, final String cigarStr, final int nm, final boolean forwardStrand)
    {
        if(!mSegmentsByAltContig.containsKey(contig))
        {
            // an unknown alt contig must not pass through as ref - that would leak _tx contig names into the BAM.
            if(contig.endsWith(ALT_CONTIG_SUFFIX))
            {
                throw new IllegalStateException(String.format("alt contig(%s) missing from the sidecar", contig));
            }

            return new LiftedAlignment(contig, pos, cigarStr, nm, false, forwardStrand, null);
        }

        ContigEntry entry = findSegment(contig, pos);
        ContigTranslateResult translated = translate(entry, pos, CigarUtils.cigarFromStr(cigarStr));
        if(translated == null)
        {
            return null;
        }

        String genomicCigar = new Cigar(mergeDeletionsIntoSplice(translated.genomicCigar().getCigarElements())).toString();
        return new LiftedAlignment(
                translated.chromosome(), translated.genomicStart(), genomicCigar, nm,
                true, forwardStrand, Orientation.fromByte((byte) entry.strand()));
    }

    public List<LiftedAlignment> liftXaAlignments(final String xaTag)
    {
        List<LiftedAlignment> lifted = new ArrayList<>();
        if(xaTag == null)
        {
            return lifted;
        }

        Set<AlignmentKey> seenKeys = new HashSet<>();
        for(String entry : xaTag.split(";"))
        {
            String[] fields = entry.split(",");
            if(fields.length != 4)
            {
                throw new IllegalStateException(String.format("malformed XA(%s) entry(%s)", xaTag, entry));
            }

            int signedPosition = Integer.parseInt(fields[1]);
            LiftedAlignment alignment = liftAlignment(
                    fields[0], Math.abs(signedPosition), fields[2], Integer.parseInt(fields[3]), signedPosition >= 0);
            if(alignment != null && seenKeys.add(alignment.key()))
            {
                lifted.add(alignment);
            }
        }

        return lifted;
    }

    // Segment owning altPos, or null when the contig is unknown. A spacer or leading-overhang position resolves to the
    // nearest segment so the overhang clamp can salvage it.
    ContigEntry findSegment(final String altContig, final int altPos)
    {
        List<ContigEntry> segments = mSegmentsByAltContig.get(altContig);
        int placement = floorIndexByContigStart(segments, altPos);
        ContigEntry segment = segments.get(placement);
        if(altPos <= segment.contigRegion().end())
        {
            return segment;
        }

        if(placement + 1 < segments.size())
        {
            ContigEntry next = segments.get(placement + 1);
            int leadingOverhang = next.contigRegion().start() - altPos;
            int trailingOverhang = altPos - segment.contigRegion().end();
            return leadingOverhang <= trailingOverhang ? next : segment;
        }

        return segment;
    }

    // index of the last segment with contigStart <= altPos, or 0 when altPos precedes the first. Mirrors
    // BaseRegion.binarySearch, which cannot be applied directly since segments are ContigEntry.
    private static int floorIndexByContigStart(final List<ContigEntry> segments, final int altPos)
    {
        int low = 0;
        int high = segments.size() - 1;
        int floor = 0;
        while(low <= high)
        {
            int mid = (low + high) >>> 1;
            if(segments.get(mid).contigRegion().start() <= altPos)
            {
                floor = mid;
                low = mid + 1;
            }
            else
            {
                high = mid - 1;
            }
        }
        return floor;
    }

    public static ContigTranslateResult translate(final ContigEntry contig, final int contigPos, final Cigar contigCigar)
    {
        List<BaseRegion> spans = contig.exonSpans();
        int clampedContigPos = contigPos;
        Cigar clampedCigar = contigCigar;

        if(contigPos < contig.contigRegion().start())
        {
            clampedCigar = clampLeadingReferenceToSoftClip(contigCigar, contig.contigRegion().start() - contigPos);
            if(clampedCigar == null)
            {
                return null;
            }

            clampedContigPos = contigPos + contigCigar.getReferenceLength() - clampedCigar.getReferenceLength();
        }

        int readEnd = clampedContigPos + clampedCigar.getReferenceLength() - 1;
        if(readEnd > contig.contigRegion().end())
        {
            clampedCigar = clampTrailingReferenceToSoftClip(clampedCigar, readEnd - contig.contigRegion().end());
            if(clampedCigar == null)
            {
                return null;
            }
        }

        int localPos = clampedContigPos - contig.contigRegion().start() + 1;
        SpanLocation start = locateStartSpan(spans, localPos);
        return walkCigarToGenome(contig, spans, start, clampedCigar);
    }

    // localPos is a 1-based offset into the concatenated exons.
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
        throw new IllegalStateException(String.format("localPos(%d) beyond exon length(%d)", localPos, exonLengthSoFar));
    }

    // Walks the clamped contig-space CIGAR, emitting an N at every exon boundary crossed.
    private static ContigTranslateResult walkCigarToGenome(
            final ContigEntry contig, final List<BaseRegion> spans, final SpanLocation start, final Cigar cigar)
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
                new Cigar(normalize(outElements)));
    }

    // A D straddling an exon boundary lifts as xD nN yD; folding the small flanking Ds into the N preserves both spans.
    static List<CigarElement> mergeDeletionsIntoSplice(final List<CigarElement> elements)
    {
        if(elements.size() < 3)
        {
            return elements;
        }

        List<CigarElement> result = new ArrayList<>(elements.size());
        for(int i = 0; i < elements.size(); ++i)
        {
            CigarElement element = elements.get(i);

            if(element.getOperator() != CigarOperator.N)
            {
                result.add(element);
                continue;
            }

            int splicedLength = element.getLength();

            // Absorb only while an aligned block still flanks that side of the splice. Absorbing a terminal deletion
            // would leave the N at the end of the cigar, which downstream reads as a junction with no exon beside it.
            if(result.size() >= 2 && isAbsorbableDeletion(result.get(result.size() - 1))
                    && isMatchOrEqualOp(result.get(result.size() - 2).getOperator()))
            {
                splicedLength += result.remove(result.size() - 1).getLength();
            }

            int afterDeletions = i + 1;
            while(afterDeletions < elements.size() && isAbsorbableDeletion(elements.get(afterDeletions)))
            {
                ++afterDeletions;
            }

            if(afterDeletions > i + 1 && afterDeletions < elements.size()
                    && isMatchOrEqualOp(elements.get(afterDeletions).getOperator()))
            {
                while(i + 1 < afterDeletions)
                {
                    splicedLength += elements.get(++i).getLength();
                }
            }

            result.add(new CigarElement(splicedLength, CigarOperator.N));
        }
        return result;
    }

    private static boolean isAbsorbableDeletion(final CigarElement element)
    {
        return element.getOperator() == CigarOperator.D && element.getLength() <= MAX_MERGED_DELETION_BP;
    }

    private record SpanLocation(int spanIndex, int genomicPos)
    {
    }

}
