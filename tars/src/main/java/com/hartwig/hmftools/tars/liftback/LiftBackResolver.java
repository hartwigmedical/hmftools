package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.common.TarsConstants.ALT_CONTIG_SUFFIX;
import static com.hartwig.hmftools.tars.common.TarsConstants.ANNOTATED_JUNCTION_MIN_ANCHOR_BP;
import static com.hartwig.hmftools.tars.common.TarsConstants.ANNOTATED_JUNCTION_MIN_SOFTCLIP_ANCHOR_BP;
import static com.hartwig.hmftools.tars.common.TarsConstants.CONFIDENT_MAPQ;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.XA_ATTRIBUTE;
import static com.hartwig.hmftools.tars.common.TarsConfig.TARS_LOGGER;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.tars.common.ContigEntry;
import com.hartwig.hmftools.tars.common.TarsConstants;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;

// Lifts SAMRecords aligned against ref + Tx contigs to genomic coords. Every input record produces exactly one
// result; UNMAPPED / LIFT_FAILED results carry placeholder fields so the emit step can flag records unmapped.
public class LiftBackResolver
{
    private static final String AS_TAG = "AS";
    private static final String XS_TAG = "XS";
    private static final String NM_TAG = "NM";

    // per-alt-contig list of segments sorted by altStart for binary search back to the owning transcript.
    private final Map<String, List<ContigEntry>> mSegmentsByAltContig;

    // when present, resolves hidden ties (XS==AS, no XA) on ref-only primaries landing inside an annotated exon.
    private final ExonRegionIndex mExonIndex;

    public LiftBackResolver(final List<ContigEntry> entries)
    {
        this(entries, null);
    }

    public LiftBackResolver(final List<ContigEntry> entries, final ExonRegionIndex exonIndex)
    {
        mSegmentsByAltContig = new HashMap<>();
        for(final ContigEntry entry : entries)
        {
            mSegmentsByAltContig.computeIfAbsent(entry.contigName(), k -> new ArrayList<>()).add(entry);
        }
        for(final List<ContigEntry> segments : mSegmentsByAltContig.values())
        {
            segments.sort(Comparator.comparingInt(ContigEntry::altStart));
        }

        mExonIndex = exonIndex;
    }

    public Set<String> contigNames()
    {
        return mSegmentsByAltContig.keySet();
    }

    // Returns the segment owning altPos, or null if it falls in an inter-transcript spacer or the contig is unknown.
    ContigEntry findSegment(final String altContig, final int altPos)
    {
        List<ContigEntry> segments = mSegmentsByAltContig.get(altContig);
        if(segments == null)
        {
            return null;
        }

        int lo = 0;
        int hi = segments.size() - 1;
        int candidate = -1;
        while(lo <= hi)
        {
            int mid = (lo + hi) >>> 1;
            if(segments.get(mid).altStart() <= altPos)
            {
                candidate = mid;
                lo = mid + 1;
            }
            else
            {
                hi = mid - 1;
            }
        }

        if(candidate < 0)
        {
            // altPos is before the first segment; return it so ContigTranslator's leading-overhang clamp can salvage it.
            return segments.isEmpty() ? null : segments.get(0);
        }

        ContigEntry segment = segments.get(candidate);
        if(altPos <= segment.altEnd())
        {
            return segment;
        }

        // altPos in the inter-segment spacer; choose the nearer neighbour and let ContigTranslator clamp the overhang.
        if(candidate + 1 < segments.size())
        {
            ContigEntry next = segments.get(candidate + 1);
            int leadingOverhang = next.altStart() - altPos;
            int trailingOverhang = altPos - segment.altEnd();
            return leadingOverhang <= trailingOverhang ? next : segment;
        }

        return segment;
    }

    // Lift-only entry point for callers that don't need the full LiftBackResult (e.g. SA tag rewriting).
    public LiftedCoords liftCoords(final String contig, final int pos, final String cigarStr)
    {
        LiftedAlignment lifted = liftAlignment(
                LiftedAlignment.AlignmentSource.SELF, contig, pos, cigarStr, 0, 0, true);
        if(lifted == null)
        {
            return null;
        }
        return new LiftedCoords(lifted.LiftedChrom, lifted.LiftedPos, lifted.LiftedCigar);
    }

    private LiftedAlignment liftSelf(final SAMRecord record)
    {
        return liftAlignment(
                LiftedAlignment.AlignmentSource.SELF,
                record.getReferenceName(), record.getAlignmentStart(), record.getCigarString(),
                getInt(record, AS_TAG), getInt(record, NM_TAG), !record.getReadNegativeStrandFlag());
    }

    // Per-call hook to normalize each candidate's lifted cigar (collapse + tail-extend) before the
    // discriminator runs, so ref-vs-tx is decided on tested facts. Supplied by the per-worker group
    // processor which owns the (thread-unsafe) engines; null on the lift-only / no-ref paths.
    @FunctionalInterface
    public interface AlignmentNormalizer
    {
        void normalize(List<LiftedAlignment> alignments, SAMRecord record);
    }

    // No-reconcile convenience for non-discriminating callers only: supplementaries and unmapped records (which
    // skip the discriminator entirely), the lift-only paths, and tests. A primary resolved through here would be
    // discriminated on raw bwa cigars, so the production primary path always uses the two-arg form with the
    // per-worker reconciler (see LiftBackGroupProcessor.decideMateGroup).
    public LiftBackResult resolve(final SAMRecord record)
    {
        return resolve(record, null);
    }

    public LiftBackResult resolve(final SAMRecord record, final AlignmentNormalizer normalizer)
    {
        if(record.getReadUnmappedFlag())
        {
            return unmappedResult(record);
        }

        if(record.getSupplementaryAlignmentFlag())
        {
            return liftSupplementary(record);
        }

        return resolvePrimary(record, normalizer);
    }

    private LiftBackResult resolvePrimary(final SAMRecord record, final AlignmentNormalizer normalizer)
    {
        List<LiftedAlignment> xaAlts = parseAndLiftXa(record);
        return resolvePrimaryWithAlts(record, xaAlts, xaAlts.size(), normalizer);
    }

    private LiftBackResult resolvePrimaryWithAlts(
            final SAMRecord record, final List<LiftedAlignment> alts, final int numXaAltsForReport,
            final AlignmentNormalizer normalizer)
    {
        LiftedAlignment self = liftSelf(record);

        if(self == null)
        {
            return unliftableResult(LiftBackResult.RecordRole.PRIMARY, numXaAltsForReport, "primary_translate_failed");
        }

        List<LiftedAlignment> allAlignments = new ArrayList<>(1 + alts.size());
        allAlignments.add(self);
        allAlignments.addAll(alts);

        // Normalize each candidate's cigar (collapse fabricated micro-junctions, reclaim terminal softclips)
        // before discriminating, so RefSoftClipped / RefFullMatch / TxHasNCigar are measured, not assumed.
        if(normalizer != null)
        {
            normalizer.normalize(allAlignments, record);
        }

        // normalize() can replace self (index 0) with a reconciled copy, so re-fetch it before marking the
        // primary choice. Marking the pre-normalize object would leave IsPrimaryChoice on an orphan no longer
        // in the list, and the discriminator's self-identity checks (apply / drop-vs-swap) would misfire.
        self = allAlignments.get(0);
        self.IsPrimaryChoice = true;

        LiftBackDiscriminator.Features features = LiftBackDiscriminator.categorize(allAlignments);
        int inputMapq = record.getMappingQuality();
        int seed = readSeed(record.getReadName());
        LiftBackDiscriminator.ApplyResult outcome = LiftBackDiscriminator.apply(
                allAlignments, features.Feature, self, seed, inputMapq != 0);
        LiftedAlignment effectivePrimary = outcome.effectivePrimary();

        List<LiftedAlignment> keptAlignments = allAlignments.stream()
                .filter(la -> !la.Dropped)
                .toList();

        int numLoci = countDistinctLoci(keptAlignments, effectivePrimary);
        int cigarsAtPrimaryLocus = countDistinctCigarsAtLocus(keptAlignments, effectivePrimary);
        String geneIds = joinGeneIds(keptAlignments);

        boolean swapped = effectivePrimary != self;
        boolean hiddenTie = inputMapq == 0 && hasHiddenTie(record);
        boolean inAnnotatedExon = mExonIndex != null
                && mExonIndex.contains(effectivePrimary.LiftedChrom, effectivePrimary.LiftedPos);
        int updatedMapq = decidePrimaryMapq(
                inputMapq, numLoci, hiddenTie, effectivePrimary.fromTxContig(), inAnnotatedExon, features.Feature);

        LiftedAlignment primaryCoords = effectivePrimary;

        if(swapped)
        {
            TARS_LOGGER.debug2("discriminator {} {}: primary -> {}:{} {} ({})",
                    features.Feature, record.getReadName(), primaryCoords.LiftedChrom, primaryCoords.LiftedPos,
                    primaryCoords.LiftedCigar, outcome.note());
        }

        return new LiftBackResult(
                RecordState.RESOLVED, features.Feature, swapped,
                LiftBackResult.Composition.fromAlignments(keptAlignments),
                LiftBackResult.RecordRole.PRIMARY,
                primaryCoords.LiftedChrom, primaryCoords.LiftedPos, primaryCoords.LiftedCigar,
                !primaryCoords.ForwardStrand,
                primaryCoords.cigarHasN(), inputMapq, updatedMapq,
                numXaAltsForReport, features.NumRefAlts, features.NumTxAlts,
                numLoci, cigarsAtPrimaryLocus,
                features.TxHasNCigar, features.TxSoftClipAtBoundary,
                features.RefSoftClipped, features.RefFullMatch,
                geneIds, outcome.note(),
                primaryCoords.TranscriptStrand,
                allAlignments);
    }

    // A supplementary is only lifted, never discriminated: lift its own coords (no XA parse, no locus pick).
    // MAPQ=0 on a tx-contig supplementary is the multi-alt-contig tie artefact; bump to CONFIDENT_MAPQ.
    private LiftBackResult liftSupplementary(final SAMRecord record)
    {
        LiftedAlignment lifted = liftSelf(record);

        if(lifted == null)
        {
            return unliftableResult(LiftBackResult.RecordRole.SUPPLEMENTARY, 0, "supp_translate_failed");
        }

        lifted.IsPrimaryChoice = true;

        int inputMapq = record.getMappingQuality();
        int outputMapq = (lifted.fromTxContig() && inputMapq == 0) ? CONFIDENT_MAPQ : inputMapq;
        int numRefAlts = lifted.fromTxContig() ? 0 : 1;
        int numTxAlts = lifted.fromTxContig() ? 1 : 0;

        return new LiftBackResult(
                RecordState.SUPPLEMENTARY, null, false,
                LiftBackResult.Composition.fromAlignments(List.of(lifted)),
                LiftBackResult.RecordRole.SUPPLEMENTARY,
                lifted.LiftedChrom, lifted.LiftedPos, lifted.LiftedCigar,
                !lifted.ForwardStrand,
                lifted.cigarHasN(), inputMapq, outputMapq,
                0, numRefAlts, numTxAlts,
                1, 1,
                false, false, false, false,
                lifted.GeneId != null ? lifted.GeneId : "", "",
                lifted.TranscriptStrand,
                List.of(lifted));
    }

    // Self is excluded from the dedup key set so a Tx XA alt lifting to the same coords as a ref self is preserved (drives BOTH_AGREE).
    private List<LiftedAlignment> parseAndLiftXa(final SAMRecord record)
    {
        List<LiftedAlignment> alts = new ArrayList<>();
        String xa = record.getStringAttribute(XA_ATTRIBUTE);
        if(xa == null || xa.isEmpty())
        {
            return alts;
        }

        Set<String> seenKeys = new HashSet<>();

        for(final String entry : xa.split(";"))
        {
            if(entry.isEmpty())
                continue;
            String[] parts = entry.split(",");
            if(parts.length < 4)
                continue;

            String contig = parts[0];
            String cigar = parts[2];

            // bwa XA pos is signed (sign encodes strand). Tolerate malformed entries.
            int signedPos;
            int nm;
            try
            {
                signedPos = Integer.parseInt(parts[1]);
            }
            catch(NumberFormatException e)
            {
                continue;
            }
            try
            {
                nm = Integer.parseInt(parts[3]);
            }
            catch(NumberFormatException e)
            {
                nm = 0;
            }

            boolean forwardStrand = signedPos >= 0;
            int pos = Math.abs(signedPos);

            LiftedAlignment lifted = liftAlignment(
                    LiftedAlignment.AlignmentSource.XA_INPUT, contig, pos, cigar, 0, nm, forwardStrand);
            if(lifted == null)
                continue;
            if(seenKeys.add(liftedKey(lifted)))
            {
                alts.add(lifted);
            }
        }

        return alts;
    }

    // Returns null when the contig is an unknown alt contig or the position falls outside any transcript.
    // Ref alignments pass through as-is.
    private LiftedAlignment liftAlignment(
            final LiftedAlignment.AlignmentSource source, final String contig, final int pos, final String cigarStr,
            final int as, final int nm, final boolean forwardStrand)
    {
        if(!mSegmentsByAltContig.containsKey(contig))
        {
            // Unknown alt contig must not pass through as ref - that would leak _tx contig names into the BAM.
            if(contig.endsWith(ALT_CONTIG_SUFFIX))
            {
                return null;
            }

            return new LiftedAlignment(
                    source, contig, pos, cigarStr,
                    contig, pos, cigarStr,
                    as, nm,
                    null, null, null,
                    false, forwardStrand);
        }

        ContigEntry entry = findSegment(contig, pos);
        if(entry == null)
        {
            return null;
        }

        Cigar parsedCigar = TextCigarCodec.decode(cigarStr);
        ContigTranslator.TranslationResult translated = ContigTranslator.translate(entry, pos, parsedCigar);
        if(translated == null)
        {
            return null;
        }

        boolean softClipAtBoundary = ContigTranslator.hasSoftClipAtExonBoundary(entry, pos, parsedCigar);

        // Drop sub-threshold terminal anchors that bwa fabricated by walking past an exon boundary.
        ContigTranslator.MicroAnchorResult trimmed = ContigTranslator.trimMicroAnchors(
                translated.genomicCigar(), ANNOTATED_JUNCTION_MIN_ANCHOR_BP, ANNOTATED_JUNCTION_MIN_SOFTCLIP_ANCHOR_BP);

        return new LiftedAlignment(
                source, contig, pos, cigarStr,
                translated.chromosome(), translated.genomicStart() + trimmed.StartShift,
                trimmed.AdjustedCigar.toString(),
                as, nm,
                entry.transName(), entry.geneId(), entry.geneName(),
                softClipAtBoundary, forwardStrand, entry.strand());
    }

    // Counts distinct genomic loci among the kept alignments, using the SAME same-locus test as buildLiftedXa
    // (overlapsPrimary) so NH and the emitted XA always agree. The primary is one locus; an alt whose genomic span
    // OVERLAPS the primary's is the same locus and collapses -- a 5'/3'-softclipped isoform copy begins at a
    // downstream exon of the same placement, lifting to a different start but a span nested inside the primary's, so
    // keying on exact start would over-count it and wrongly withhold the MAPQ bump. Alts that do NOT overlap the
    // primary are genuinely distinct placements (repeats/paralogs, other chromosomes); they are interval-merged
    // among themselves so a cluster of overlapping repeat copies at one other locus counts once, but a tandem-repeat
    // array offset beyond the primary is NOT chained back into the primary (that would falsely bump a repeat
    // multimapper). numLoci therefore equals 1 exactly when every alt overlaps the primary -- i.e. when buildLiftedXa
    // emits no XA -- so a bumped MAPQ-60 / NH-1 read never carries a stale XA.
    private static int countDistinctLoci(final List<LiftedAlignment> alignments, final LiftedAlignment primary)
    {
        int primaryEnd = spanEnd(primary);

        Map<String, List<int[]>> distinctSpans = new HashMap<>();
        distinctSpans.computeIfAbsent(primary.LiftedChrom, k -> new ArrayList<>())
                .add(new int[] { primary.LiftedPos, primaryEnd });

        for(final LiftedAlignment la : alignments)
        {
            if(la == primary)
            {
                continue;
            }
            int end = spanEnd(la);
            boolean overlapsPrimary = la.LiftedChrom.equals(primary.LiftedChrom)
                    && la.LiftedPos <= primaryEnd && primary.LiftedPos <= end;
            if(overlapsPrimary)
            {
                continue;
            }
            distinctSpans.computeIfAbsent(la.LiftedChrom, k -> new ArrayList<>()).add(new int[] { la.LiftedPos, end });
        }

        int loci = 0;
        for(final List<int[]> spans : distinctSpans.values())
        {
            spans.sort(Comparator.comparingInt(s -> s[0]));
            int clusterEnd = -1;
            for(final int[] span : spans)
            {
                if(span[0] > clusterEnd)
                {
                    ++loci;
                    clusterEnd = span[1];
                }
                else
                {
                    clusterEnd = Math.max(clusterEnd, span[1]);
                }
            }
        }
        return loci;
    }

    private static int spanEnd(final LiftedAlignment la)
    {
        return la.LiftedPos + TextCigarCodec.decode(la.LiftedCigar).getReferenceLength() - 1;
    }

    private static int countDistinctCigarsAtLocus(final List<LiftedAlignment> alignments, final LiftedAlignment primary)
    {
        String primaryLocus = locusKey(primary);
        Set<String> cigars = new HashSet<>();
        for(final LiftedAlignment la : alignments)
        {
            if(locusKey(la).equals(primaryLocus))
            {
                cigars.add(la.LiftedCigar);
            }
        }
        return cigars.size();
    }

    private static String joinGeneIds(final List<LiftedAlignment> alignments)
    {
        Set<String> geneIds = new HashSet<>();
        for(final LiftedAlignment la : alignments)
            if(la.GeneId != null)
            {
                geneIds.add(la.GeneId);
            }
        return geneIds.stream().sorted().collect(Collectors.joining("|"));
    }

    private static String liftedKey(final LiftedAlignment la)
    {
        return la.LiftedChrom + ":" + la.LiftedPos + ":" + la.LiftedCigar;
    }

    private static String locusKey(final LiftedAlignment la)
    {
        return la.LiftedChrom + ":" + la.LiftedPos;
    }

    private static int getInt(final SAMRecord record, final String tag)
    {
        Integer val = record.getIntegerAttribute(tag);
        return val != null ? val : 0;
    }

    private static LiftBackResult unmappedResult(final SAMRecord record)
    {
        LiftBackResult.RecordRole role = record.isSecondaryOrSupplementary()
                ? LiftBackResult.RecordRole.SUPPLEMENTARY
                : LiftBackResult.RecordRole.PRIMARY;

        return new LiftBackResult(
                RecordState.UNMAPPED, null, false, LiftBackResult.Composition.NONE,
                role,
                "*", 0, "*",
                false,
                false, 0, 0,
                0, 0, 0,
                0, 0,
                false, false, false, false,
                "", "",
                0,
                List.of());
    }

    // Bump a MAPQ-0 primary to 60 only when it lifts to a single locus with no competing alt. Multi-locus reads
    // (including a discriminator swap to an alt at another locus) are never promoted. A hidden tie (XS==AS, an
    // equal-scoring alt bwa did not emit) is a latent second locus and blocks the bump unless tx provenance or an
    // annotated exon vouches for the placement.
    static int decidePrimaryMapq(
            final int inputMapq, final int numLoci, final boolean hiddenTie,
            final boolean primaryFromTxContig, final boolean primaryInAnnotatedExon)
    {
        return decidePrimaryMapq(inputMapq, numLoci, hiddenTie, primaryFromTxContig, primaryInAnnotatedExon, null);
    }

    static int decidePrimaryMapq(
            final int inputMapq, final int numLoci, final boolean hiddenTie,
            final boolean primaryFromTxContig, final boolean primaryInAnnotatedExon, final DecidingFeature feature)
    {
        if(inputMapq != 0)
        {
            return inputMapq;
        }
        // an ambiguous single-locus call is a coin-flip between tx and ref, so it must never be promoted to confident.
        if(feature == DecidingFeature.AMBIGUOUS)
        {
            return inputMapq;
        }
        if(numLoci != 1)
        {
            return inputMapq;
        }
        boolean unresolvedHiddenTie = hiddenTie && !primaryFromTxContig && !primaryInAnnotatedExon;
        if(unresolvedHiddenTie)
        {
            return inputMapq;
        }
        return CONFIDENT_MAPQ;
    }

    // Per-read deterministic seed for the not-confident random placement picks. Both mates share a read name,
    // so a pair is placed together, and the value is stable across runs, keeping the output reproducible.
    static int readSeed(final String readName)
    {
        int hash = readName.hashCode();
        return hash ^ (hash >>> 16);
    }

    // When XS == AS, an equally-scoring alt wasn't emitted by bwa. Flag as a hidden tie to skip the MAPQ bump.
    private static boolean hasHiddenTie(final SAMRecord record)
    {
        Integer alignmentScore = record.getIntegerAttribute(AS_TAG);
        Integer suboptimalScore = record.getIntegerAttribute(XS_TAG);
        return alignmentScore != null && suboptimalScore != null && suboptimalScore.intValue() == alignmentScore.intValue();
    }

    private static LiftBackResult unliftableResult(final LiftBackResult.RecordRole role, final int numXaAlts, final String note)
    {
        return new LiftBackResult(
                RecordState.LIFT_FAILED, null, false, LiftBackResult.Composition.NONE,
                role,
                "*", 0, "*",
                false,
                false, 0, 0,
                numXaAlts, 0, 0,
                0, 0,
                false, false, false, false,
                "", note,
                0,
                List.of());
    }
}
