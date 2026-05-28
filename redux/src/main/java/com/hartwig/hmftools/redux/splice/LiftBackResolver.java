package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.redux.splice.SpliceCommon.ALT_CONTIG_SUFFIX;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;

// per-record categorizer: takes a SAMRecord aligned against ref + Tx contigs and produces a LiftBackResult
// with the chosen lifted (chrom, pos, CIGAR), the assigned LiftBackCategory, and the full lifted alignment set.
//
// 1:1 contract: every input record produces exactly one result. UNMAPPED / LIFT_FAILED results carry placeholder
// fields so the downstream emit step can flag the BAM record unmapped without dropping it.
public class LiftBackResolver
{
    private static final String XA_TAG = "XA";
    private static final String AS_TAG = "AS";
    private static final String XS_TAG = "XS";
    private static final String NM_TAG = "NM";

    private static final int RESCUED_MAPQ = 60;
    private static final int INPUT_UNIQUE_MAPQ = 60;

    // Tail-anchor floor for tx-contig supps. A trailing yM after an N with y below this is dropped
    // into the trailing softclip — the junction it implies isn't supported by enough matched bases.
    private static final int ANNOTATED_JUNCTION_MIN_ANCHOR_BP = 3;

    // per-alt-contig list of segments sorted by altStart so a record's alt-contig position can be
    // bin-searched back to the owning transcript. Non-alt-contig alignments (ref) fall through unindexed.
    private final Map<String, List<ContigEntry>> mSegmentsByAltContig;

    // when true, the input BAM was produced with bwa-mem2 -a so alts arrive as secondary records (0x100)
    // instead of via the XA tag. Skips XA parsing on primaries and disables the tx-contig MAPQ rescue on
    // secondaries (their MAPQ=0 is a real ambiguity signal, not the multi-alt-contig tie artefact).
    private final boolean mEmitSecondaries;

    // optional annotated-exon lookup. When provided, a hidden tie (XS==AS, no XA) on a ref-only primary
    // is overridden when the lifted primary lands inside an annotated exon — the tied alt is then
    // almost certainly a sub-threshold intronic/intergenic paralog.
    private final ExonRegionIndex mExonIndex;

    public LiftBackResolver(final List<ContigEntry> entries)
    {
        this(entries, false, null);
    }

    public LiftBackResolver(final List<ContigEntry> entries, final boolean emitSecondaries)
    {
        this(entries, emitSecondaries, null);
    }

    public LiftBackResolver(final List<ContigEntry> entries, final boolean emitSecondaries, final ExonRegionIndex exonIndex)
    {
        mSegmentsByAltContig = new HashMap<>();
        for(final ContigEntry entry : entries)
            mSegmentsByAltContig.computeIfAbsent(entry.contigName(), k -> new ArrayList<>()).add(entry);
        for(final List<ContigEntry> segments : mSegmentsByAltContig.values())
            segments.sort(Comparator.comparingInt(ContigEntry::altStart));

        mEmitSecondaries = emitSecondaries;
        mExonIndex = exonIndex;
    }

    public Set<String> contigNames()
    {
        return mSegmentsByAltContig.keySet();
    }

    // locates the transcript segment that owns altPos on the given alt contig, or null if the position
    // falls outside any transcript (e.g. inside the inter-transcript N spacer) or the contig isn't an alt contig.
    ContigEntry findSegment(final String altContig, final int altPos)
    {
        final List<ContigEntry> segments = mSegmentsByAltContig.get(altContig);
        if(segments == null)
            return null;

        int lo = 0;
        int hi = segments.size() - 1;
        int candidate = -1;
        while(lo <= hi)
        {
            final int mid = (lo + hi) >>> 1;
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
            // altPos sits before the first segment's altStart (in the upstream spacer). Hand back the first
            // segment so ContigTranslator's leading-overhang clamp can salvage it.
            return segments.isEmpty() ? null : segments.get(0);
        }

        final ContigEntry segment = segments.get(candidate);
        if(altPos <= segment.altEnd())
            return segment;

        // altPos sits in the spacer between candidate and candidate+1. Choose whichever neighbour the read
        // overhangs less, and let ContigTranslator's clamp convert the overhang into soft-clip.
        if(candidate + 1 < segments.size())
        {
            final ContigEntry next = segments.get(candidate + 1);
            final int leadingOverhang = next.altStart() - altPos;
            final int trailingOverhang = altPos - segment.altEnd();
            return leadingOverhang <= trailingOverhang ? next : segment;
        }

        return segment;
    }

    // lift-only API for callers that don't need the full LiftBackResult machinery (e.g. SA tag rewriting).
    public LiftedCoords liftCoords(final String contig, final int pos, final String cigarStr)
    {
        final LiftedAlignment lifted = liftAlignment(
                LiftedAlignment.AlignmentSource.SELF, contig, pos, cigarStr, 0, 0, true);
        if(lifted == null)
            return null;
        return new LiftedCoords(lifted.LiftedChrom, lifted.LiftedPos, lifted.LiftedCigar);
    }

    public LiftBackResult resolve(final SAMRecord record)
    {
        if(record.getReadUnmappedFlag())
            return unmappedResult(record);

        if(record.getSupplementaryAlignmentFlag())
            return supplementaryResult(record);

        if(record.isSecondaryAlignment())
            return nonPrimaryResult(record, NonPrimaryKind.SECONDARY);

        return resolvePrimary(record);
    }

    // Non-primary records (secondaries from -emit_secondaries and supplementaries from split reads)
    // share their lift/result shape; only the MAPQ-rescue rule and failure-note string differ.
    enum NonPrimaryKind
    {
        // 0x100. MAPQ=0 means real ambiguity at this locus; never rescue.
        SECONDARY("secondary_translate_failed", false),
        // 0x800. MAPQ=0 on a tx-contig supp is the multi-alt-contig tie artefact; rescue it.
        SUPPLEMENTARY("supp_translate_failed", true);

        final String FailureNote;
        final boolean RescueZeroMapqOnTxContig;

        NonPrimaryKind(final String failureNote, final boolean rescueZeroMapqOnTxContig)
        {
            FailureNote = failureNote;
            RescueZeroMapqOnTxContig = rescueZeroMapqOnTxContig;
        }
    }

    private LiftBackResult resolvePrimary(final SAMRecord record)
    {
        // in -emit_secondaries mode the alts arrive as separate 0x100 records, so any XA entries on the
        // primary (if bwa-mem2 still emits some) would double-count. Skip the parse — the per-pair grouped
        // path (see resolvePrimaryWithSecondaries) is what feeds the real alt set in that mode.
        final List<LiftedAlignment> xaAlts = mEmitSecondaries ? List.of() : parseAndLiftXa(record);
        // numXaAlts in the result tracks the *deduped + lifted* count, not the raw XA entry count, so dedup
        // behavior remains visible in TSV-A.
        return resolvePrimaryWithAlts(record, xaAlts, xaAlts.size());
    }

    // grouped-mode entry point: lifts each secondary as a SECONDARY_INPUT alt and feeds them into the
    // primary's alignment set. Used when bwa-mem2 -a put alts on separate 0x100 records but the BAM has
    // been name-sorted so a primary + its secondaries can be processed as one group.
    public LiftBackResult resolvePrimaryWithSecondaries(final SAMRecord primary, final List<SAMRecord> secondaryRecords)
    {
        final List<LiftedAlignment> alts = new ArrayList<>(secondaryRecords.size());
        final Set<String> seenKeys = new HashSet<>();
        for(final SAMRecord secondary : secondaryRecords)
        {
            if(secondary.getReadUnmappedFlag())
                continue;
            final LiftedAlignment lifted = liftAlignment(
                    LiftedAlignment.AlignmentSource.SECONDARY_INPUT,
                    secondary.getReferenceName(), secondary.getAlignmentStart(), secondary.getCigarString(),
                    getInt(secondary, AS_TAG), getInt(secondary, NM_TAG),
                    !secondary.getReadNegativeStrandFlag());
            if(lifted == null)
                continue;
            if(seenKeys.add(liftedKey(lifted)))
                alts.add(lifted);
        }
        return resolvePrimaryWithAlts(primary, alts, secondaryRecords.size());
    }

    // shared body for both XA-mode (alts from primary's XA tag) and secondaries-mode (alts from sibling
    // 0x100 records). numXaAltsForReport is what the LiftBackResult.numXaAlts field carries downstream —
    // XA-mode passes the XA tag count, secondaries-mode passes the secondary-record count.
    private LiftBackResult resolvePrimaryWithAlts(
            final SAMRecord record, final List<LiftedAlignment> alts, final int numXaAltsForReport)
    {
        final LiftedAlignment self = liftAlignment(
                LiftedAlignment.AlignmentSource.SELF,
                record.getReferenceName(), record.getAlignmentStart(), record.getCigarString(),
                getInt(record, AS_TAG), getInt(record, NM_TAG),
                !record.getReadNegativeStrandFlag());

        if(self == null)
            return unliftableResult(LiftBackResult.RecordRole.PRIMARY, numXaAltsForReport, "primary_translate_failed");

        self.IsPrimaryChoice = true;

        final List<LiftedAlignment> allAlignments = new ArrayList<>(1 + alts.size());
        allAlignments.add(self);
        allAlignments.addAll(alts);

        final LiftBackDiscriminator.Features features = LiftBackDiscriminator.categorize(allAlignments);
        final LiftBackDiscriminator.Outcome outcome = LiftBackDiscriminator.apply(allAlignments, features.Category, self);
        final LiftedAlignment effectivePrimary = outcome.effectivePrimary();

        final List<LiftedAlignment> keptAlignments = allAlignments.stream()
                .filter(la -> !la.Dropped)
                .collect(Collectors.toList());

        final int numLoci = countDistinctLoci(keptAlignments);
        final int cigarsAtPrimaryLocus = countDistinctCigarsAtLocus(keptAlignments, effectivePrimary);
        final String geneIds = joinGeneIds(keptAlignments);

        final int inputMapq = record.getMappingQuality();
        final boolean swapped = effectivePrimary != self;
        final boolean hiddenTie = inputMapq == 0 && hasHiddenTie(record);
        final boolean inAnnotatedExon = mExonIndex != null
                && mExonIndex.contains(effectivePrimary.LiftedChrom, effectivePrimary.LiftedPos);
        final int updatedMapq = decidePrimaryMapq(
                inputMapq, numLoci, swapped, hiddenTie, effectivePrimary.fromTxContig(), inAnnotatedExon);

        // in secondaries mode each alt is its own SAM record and will be written separately with its own
        // lifted coords. The primary record must keep self's lifted coords (not the discriminator's
        // effective-primary swap target), otherwise the primary and the chosen-alt secondary would both
        // land at the same locus and the original primary's locus would be lost.
        final LiftedAlignment primaryCoords = mEmitSecondaries ? self : effectivePrimary;

        return new LiftBackResult(
                features.Category, LiftBackResult.Composition.fromAlignments(keptAlignments),
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

    private LiftBackResult supplementaryResult(final SAMRecord record)
    {
        return nonPrimaryResult(record, NonPrimaryKind.SUPPLEMENTARY);
    }

    // Shared body for secondaries (0x100, from -emit_secondaries) and supplementaries (0x800, from
    // split reads). Lift coords, build a SUPPLEMENTARY-role result with one alignment in the set.
    // MAPQ rescue policy is kind-dependent: see NonPrimaryKind doc.
    LiftBackResult nonPrimaryResult(final SAMRecord record, final NonPrimaryKind kind)
    {
        final LiftedAlignment lifted = liftAlignment(
                LiftedAlignment.AlignmentSource.SELF,
                record.getReferenceName(), record.getAlignmentStart(), record.getCigarString(),
                getInt(record, AS_TAG), getInt(record, NM_TAG),
                !record.getReadNegativeStrandFlag());

        if(lifted == null)
            return unliftableResult(LiftBackResult.RecordRole.SUPPLEMENTARY, 0, kind.FailureNote);

        lifted.IsPrimaryChoice = true;

        final int inputMapq = record.getMappingQuality();
        final int outputMapq = (kind.RescueZeroMapqOnTxContig && lifted.fromTxContig() && inputMapq == 0)
                ? RESCUED_MAPQ : inputMapq;
        final int numRefAlts = lifted.fromTxContig() ? 0 : 1;
        final int numTxAlts = lifted.fromTxContig() ? 1 : 0;

        return new LiftBackResult(
                LiftBackCategory.SUPPLEMENTARY, LiftBackResult.Composition.fromAlignments(List.of(lifted)),
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

    // parses the XA tag, lifts each alt, and dedupes XA-internally by lifted (chrom, pos, CIGAR).
    // Self is intentionally NOT in the dedup key set so a Tx XA alt that lifts to the same coords as a
    // ref self-alignment is preserved (drives BOTH_AGREE).
    private List<LiftedAlignment> parseAndLiftXa(final SAMRecord record)
    {
        final List<LiftedAlignment> alts = new ArrayList<>();
        final String xa = record.getStringAttribute(XA_TAG);
        if(xa == null || xa.isEmpty())
            return alts;

        final Set<String> seenKeys = new HashSet<>();

        for(final String entry : xa.split(";"))
        {
            if(entry.isEmpty())
                continue;
            final String[] parts = entry.split(",");
            if(parts.length < 4)
                continue;

            final String contig = parts[0];
            final String cigar = parts[2];

            // bwa XA pos is signed (sign = strand). Tolerate malformed XA — skip entries we can't parse.
            final int signedPos;
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

            final boolean forwardStrand = signedPos >= 0;
            final int pos = Math.abs(signedPos);

            final LiftedAlignment lifted = liftAlignment(
                    LiftedAlignment.AlignmentSource.XA_INPUT, contig, pos, cigar, 0, nm, forwardStrand);
            if(lifted == null)
                continue;
            if(seenKeys.add(liftedKey(lifted)))
                alts.add(lifted);
        }

        return alts;
    }

    // single lift kernel. Returns null when the contig is an alt contig but the position cannot be
    // translated (alt missing from segments map, or position falls outside any transcript span).
    // For ref alignments the returned LiftedAlignment is a coord-preserving pass-through.
    private LiftedAlignment liftAlignment(
            final LiftedAlignment.AlignmentSource source, final String contig, final int pos, final String cigarStr,
            final int as, final int nm, final boolean forwardStrand)
    {
        if(!mSegmentsByAltContig.containsKey(contig))
        {
            // alt contig missing from the segments map (FASTA/sidecar mismatch) must not pass through
            // as if it were ref — the resulting "lifted" coords would leak _tx contig names into the BAM.
            if(contig.endsWith(ALT_CONTIG_SUFFIX))
                return null;

            return new LiftedAlignment(
                    source, contig, pos, cigarStr,
                    contig, pos, cigarStr,
                    as, nm,
                    null, null, null,
                    false, forwardStrand);
        }

        final ContigEntry entry = findSegment(contig, pos);
        if(entry == null)
            return null;

        final Cigar parsedCigar = TextCigarCodec.decode(cigarStr);
        final ContigTranslator.TranslationResult translated = ContigTranslator.translate(entry, pos, parsedCigar);
        if(translated == null)
            return null;

        final boolean softClipAtBoundary = ContigTranslator.hasSoftClipAtExonBoundary(entry, pos, parsedCigar);

        // A tx-contig walk that dribbles a few bases past an exon boundary fabricates a junction
        // anchoring a tiny terminal yM. Drop sub-threshold anchors at both ends so we don't emit a
        // junction we can't support; a leading trim advances the start past the dropped anchor.
        final ContigTranslator.MicroAnchorResult trimmed = ContigTranslator.trimMicroAnchors(
                translated.genomicCigar(), ANNOTATED_JUNCTION_MIN_ANCHOR_BP);

        return new LiftedAlignment(
                source, contig, pos, cigarStr,
                translated.chromosome(), translated.genomicStart() + trimmed.StartShift,
                trimmed.AdjustedCigar.toString(),
                as, nm,
                entry.transName(), entry.geneId(), entry.geneName(),
                softClipAtBoundary, forwardStrand, entry.strand());
    }

    private static int countDistinctLoci(final List<LiftedAlignment> alignments)
    {
        final Set<String> loci = new HashSet<>();
        for(final LiftedAlignment la : alignments)
            loci.add(locusKey(la));
        return loci.size();
    }

    private static int countDistinctCigarsAtLocus(final List<LiftedAlignment> alignments, final LiftedAlignment primary)
    {
        final String primaryLocus = locusKey(primary);
        final Set<String> cigars = new HashSet<>();
        for(final LiftedAlignment la : alignments)
        {
            if(locusKey(la).equals(primaryLocus))
                cigars.add(la.LiftedCigar);
        }
        return cigars.size();
    }

    private static String joinGeneIds(final List<LiftedAlignment> alignments)
    {
        final Set<String> geneIds = new HashSet<>();
        for(final LiftedAlignment la : alignments)
            if(la.GeneId != null)
                geneIds.add(la.GeneId);
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
        final Integer val = record.getIntegerAttribute(tag);
        return val != null ? val : 0;
    }

    private static LiftBackResult unmappedResult(final SAMRecord record)
    {
        final LiftBackResult.RecordRole role = record.isSecondaryOrSupplementary()
                ? LiftBackResult.RecordRole.SUPPLEMENTARY
                : LiftBackResult.RecordRole.PRIMARY;

        return new LiftBackResult(
                LiftBackCategory.UNMAPPED, LiftBackResult.Composition.NONE,
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

    // MAPQ rewrite policy for primary records. Pure function of the inputs; extracted so the policy
    // can be unit-tested independently of LiftBackDiscriminator / SAMRecord plumbing.
    //
    //  (1) primary was swapped by the discriminator → rescue to RESCUED_MAPQ.
    //  (2) input MAPQ=0 and the read lifts back to a single genomic locus → rescue, unless an
    //      unresolved hidden tie is in play (XS==AS, primary isn't tx-derived, and we have no exon
    //      evidence that the unseen alt is sub-threshold).
    //  (3) input MAPQ=60 is the sanger cap for confident-unique placements → pass through as
    //      RESCUED_MAPQ (no-op when RESCUED_MAPQ == 60).
    //  (4) input MAPQ in (0, 60) is graded quality signal → leave alone.
    static int decidePrimaryMapq(
            final int inputMapq, final int numLoci, final boolean swapped, final boolean hiddenTie,
            final boolean primaryFromTxContig, final boolean primaryInAnnotatedExon)
    {
        if(swapped)
            return RESCUED_MAPQ;
        final boolean unresolvedHiddenTie = hiddenTie && !primaryFromTxContig && !primaryInAnnotatedExon;
        if(numLoci == 1 && inputMapq == 0 && !unresolvedHiddenTie)
            return RESCUED_MAPQ;
        if(inputMapq == INPUT_UNIQUE_MAPQ)
            return RESCUED_MAPQ;
        return inputMapq;
    }

    // bwa-mem2 reports XS as the second-best alignment score; when XS == AS, an equally-scoring alt
    // exists that wasn't emitted (XA omitted because of suppression heuristics, contig type, etc.). The
    // resolver can't see this alt so it would mistake the read for unique — flag it as a hidden tie so
    // the MAPQ rescue is skipped.
    private static boolean hasHiddenTie(final SAMRecord record)
    {
        final Integer alignmentScore = record.getIntegerAttribute(AS_TAG);
        final Integer suboptimalScore = record.getIntegerAttribute(XS_TAG);
        return alignmentScore != null && suboptimalScore != null && suboptimalScore.intValue() == alignmentScore.intValue();
    }

    private static LiftBackResult unliftableResult(final LiftBackResult.RecordRole role, final int numXaAlts, final String note)
    {
        return new LiftBackResult(
                LiftBackCategory.LIFT_FAILED, LiftBackResult.Composition.NONE,
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
