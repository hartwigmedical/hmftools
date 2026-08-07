package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipped;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightSoftClipped;
import static com.hartwig.hmftools.tars.common.TarsCigarUtils.clampLeadingReferenceToSoftClip;
import static com.hartwig.hmftools.tars.common.TarsCigarUtils.clampTrailingReferenceToSoftClip;
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
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.tars.common.ContigEntry;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

// Translates a transcript-contig alignment (pos + cigar) to genome coordinates, re-inserting introns as N gaps.
// Holds the sidecar's segments per alt contig so a position can be resolved to its owning transcript first;
// translate() is the coordinate walk on its own and needs no index.
public final class ContigTranslator
{
    // per-alt-contig segments sorted by contigStart for binary search back to the owning transcript.
    private final Map<String, List<ContigEntry>> mSegmentsByAltContig;

    public ContigTranslator(final List<ContigEntry> entries)
    {
        mSegmentsByAltContig = new HashMap<>();
        for(ContigEntry entry : entries)
        {
            mSegmentsByAltContig.computeIfAbsent(entry.contigName(), key -> new ArrayList<>()).add(entry);
        }
        for(List<ContigEntry> segments : mSegmentsByAltContig.values())
        {
            segments.sort(Comparator.comparingInt(ContigEntry::contigStart));
        }
    }

    public Set<String> contigNames()
    {
        return mSegmentsByAltContig.keySet();
    }

    // Lifts one placement to genomic coordinates. Returns null when the contig is an unknown alt contig or the position
    // falls outside any transcript; ref-genome placements pass through unchanged.
    public LiftedAlignment liftAlignment(
            final String contig, final int pos, final String cigarStr, final int nm, final boolean forwardStrand)
    {
        if(!mSegmentsByAltContig.containsKey(contig))
        {
            // an unknown alt contig must not pass through as ref - that would leak _tx contig names into the BAM.
            if(contig.endsWith(ALT_CONTIG_SUFFIX))
            {
                return null;
            }

            return new LiftedAlignment(contig, pos, cigarStr, nm, false, false, forwardStrand, 0);
        }

        ContigEntry entry = findSegment(contig, pos);
        if(entry == null)
        {
            return null;
        }

        ContigTranslateResult translated = translate(entry, pos, CigarUtils.cigarFromStr(cigarStr));
        if(translated == null)
        {
            return null;
        }

        String genomicCigar = new Cigar(mergeDeletionsIntoSplice(translated.genomicCigar().getCigarElements())).toString();
        return new LiftedAlignment(
                translated.chromosome(), translated.genomicStart(), genomicCigar, nm,
                true, translated.softClipAtExonBoundary(), forwardStrand, entry.strand());
    }

    // Lifts each XA entry, skips invalid placements and keeps one copy of each lifted placement.
    public List<LiftedAlignment> liftXaAlignments(final String xaTag)
    {
        List<LiftedAlignment> lifted = new ArrayList<>();
        if(xaTag == null || xaTag.isEmpty())
        {
            return lifted;
        }

        Set<String> seenKeys = new HashSet<>();
        for(String entry : xaTag.split(";"))
        {
            if(entry.isEmpty())
                continue;

            String[] fields = entry.split(",");
            if(fields.length < 4)
                continue;

            int signedPosition;
            try
            {
                signedPosition = Integer.parseInt(fields[1]);
            }
            catch(NumberFormatException e)
            {
                continue;
            }

            int numMismatches;
            try
            {
                numMismatches = Integer.parseInt(fields[3]);
            }
            catch(NumberFormatException e)
            {
                numMismatches = 0;
            }

            LiftedAlignment alignment = liftAlignment(
                    fields[0], Math.abs(signedPosition), fields[2], numMismatches, signedPosition >= 0);
            if(alignment != null && seenKeys.add(liftedKey(alignment)))
            {
                lifted.add(alignment);
            }
        }

        return lifted;
    }

    private static String liftedKey(final LiftedAlignment alignment)
    {
        return alignment.LiftedChromosome + ":" + alignment.LiftedPos + ":" + alignment.LiftedCigar
                + ":" + (alignment.ForwardStrand ? '+' : '-');
    }

    // Segment owning altPos, or null when the contig is unknown. A spacer or leading-overhang position resolves to the
    // nearest segment so the overhang clamp below can salvage it.
    ContigEntry findSegment(final String altContig, final int altPos)
    {
        List<ContigEntry> segments = mSegmentsByAltContig.get(altContig);
        if(segments == null || segments.isEmpty())
        {
            return null;
        }

        int candidate = floorIndexByContigStart(segments, altPos);
        ContigEntry segment = segments.get(candidate);
        if(altPos <= segment.contigEnd())
        {
            return segment;
        }

        // altPos in an inter-segment spacer: pick the nearer neighbour to clamp against.
        if(candidate + 1 < segments.size())
        {
            ContigEntry next = segments.get(candidate + 1);
            int leadingOverhang = next.contigStart() - altPos;
            int trailingOverhang = altPos - segment.contigEnd();
            return leadingOverhang <= trailingOverhang ? next : segment;
        }

        return segment;
    }

    // index of the last segment with contigStart <= altPos, or 0 when altPos precedes the first, mirroring
    // BaseRegion.binarySearch (segments are ContigEntry, not BaseRegion, so its overload cannot be applied directly).
    private static int floorIndexByContigStart(final List<ContigEntry> segments, final int altPos)
    {
        int low = 0;
        int high = segments.size() - 1;
        int floor = 0;
        while(low <= high)
        {
            int mid = (low + high) >>> 1;
            if(segments.get(mid).contigStart() <= altPos)
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
            cigar = clampLeadingReferenceToSoftClip(cigar, overhang);
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
            cigar = clampTrailingReferenceToSoftClip(cigar, overhang);
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
                new Cigar(normalize(outElements)),
                softClipAtExonBoundary);
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

            if(!result.isEmpty() && isAbsorbableDeletion(result.get(result.size() - 1)))
            {
                splicedLength += result.remove(result.size() - 1).getLength();
            }

            while(i + 1 < elements.size() && isAbsorbableDeletion(elements.get(i + 1)))
            {
                splicedLength += elements.get(++i).getLength();
            }

            result.add(new CigarElement(splicedLength, CigarOperator.N));
        }
        return result;
    }

    private static boolean isAbsorbableDeletion(final CigarElement element)
    {
        return element.getOperator() == CigarOperator.D && element.getLength() <= MAX_MERGED_DELETION_BP;
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
}
