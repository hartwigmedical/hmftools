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

    public LiftBackResolver(final List<ContigEntry> entries)
    {
        this(entries, false);
    }

    public LiftBackResolver(final List<ContigEntry> entries, final boolean emitSecondaries)
    {
        mSegmentsByAltContig = new HashMap<>();
        for(final ContigEntry entry : entries)
            mSegmentsByAltContig.computeIfAbsent(entry.contigName(), k -> new ArrayList<>()).add(entry);
        for(final List<ContigEntry> segments : mSegmentsByAltContig.values())
            segments.sort(Comparator.comparingInt(ContigEntry::altStart));

        mEmitSecondaries = emitSecondaries;
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
            return secondaryResult(record);

        return resolvePrimary(record);
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

        // MAPQ rewrite rules:
        //  (1) input MAPQ=0 from a ref/tx tie that lifts back to one genomic locus is an artefact of the
        //      tx-contig representation — restore to RESCUED_MAPQ so downstream tools (Isofox, redux)
        //      don't discard as ambiguous. Also rescue when we swapped the BAM primary.
        //  (2) input MAPQ=60 is the standard sanger cap for confident-unique placements; passes through
        //      unchanged (no-op now that RESCUED_MAPQ == 60). STAR's NH==1 sentinel of 255 and bwa's 60
        //      both mean "unique"; the downstream comparison treats (255, 60) as a match.
        //  (3) input MAPQ in (0, 60) carries graded alignment-quality signal (ambiguous boundaries, low
        //      complexity, indels near edges) that is strictly more informative than a binary sentinel —
        //      leave it alone.
        // Guard on (1): if the input record reports XS == AS, an equally-scoring alt is hiding (XA
        // suppressed by the source aligner, e.g. numt / chrM paralog). The chosen primary's tx provenance
        // overrides — a spliced-transcript match is stronger evidence than a sequence-only tie. Otherwise
        // keep MAPQ at 0 so downstream tools see honest ambiguity.
        // TODO extend the override: when the lifted primary falls cleanly inside an annotated exon AND
        // the hidden tie has no XA, the tie is a sub-threshold paralog and the rescue should still fire.
        // Needs an exon annotation index threaded into the resolver.
        final int inputMapq = record.getMappingQuality();
        final boolean swapped = effectivePrimary != self;
        final boolean unresolvedHiddenTie = inputMapq == 0 && hasHiddenTie(record) && !effectivePrimary.fromTxContig();
        final int updatedMapq;
        if(swapped)
            updatedMapq = RESCUED_MAPQ;
        else if(numLoci == 1 && inputMapq == 0 && !unresolvedHiddenTie)
            updatedMapq = RESCUED_MAPQ;
        else if(inputMapq == INPUT_UNIQUE_MAPQ)
            updatedMapq = RESCUED_MAPQ;
        else
            updatedMapq = inputMapq;

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
                allAlignments);
    }

    // secondary records (0x100) in -emit_secondaries mode are bwa-mem2's alt mappings expressed as their own
    // records. Lift coords like a supplementary, but DO NOT apply the tx-contig MAPQ rescue: a secondary's
    // MAPQ=0 reflects genuine ambiguity at that locus (multi-mapper), not the tied-alt-contig artefact that
    // motivates the primary's rescue. Returns a SUPPLEMENTARY-role result so the downstream writer treats it
    // the same way (no XA rewrite, mate fields patched from cache).
    private LiftBackResult secondaryResult(final SAMRecord record)
    {
        final LiftedAlignment lifted = liftAlignment(
                LiftedAlignment.AlignmentSource.SELF,
                record.getReferenceName(), record.getAlignmentStart(), record.getCigarString(),
                getInt(record, AS_TAG), getInt(record, NM_TAG),
                !record.getReadNegativeStrandFlag());

        if(lifted == null)
            return unliftableResult(LiftBackResult.RecordRole.SUPPLEMENTARY, 0, "secondary_translate_failed");

        lifted.IsPrimaryChoice = true;

        final int inputMapq = record.getMappingQuality();
        final int numRefAlts = lifted.fromTxContig() ? 0 : 1;
        final int numTxAlts = lifted.fromTxContig() ? 1 : 0;

        return new LiftBackResult(
                LiftBackCategory.SUPPLEMENTARY, LiftBackResult.Composition.fromAlignments(List.of(lifted)),
                LiftBackResult.RecordRole.SUPPLEMENTARY,
                lifted.LiftedChrom, lifted.LiftedPos, lifted.LiftedCigar,
                !lifted.ForwardStrand,
                lifted.cigarHasN(), inputMapq, inputMapq,
                0, numRefAlts, numTxAlts,
                1, 1,
                false, false, false, false,
                lifted.GeneId != null ? lifted.GeneId : "", "",
                List.of(lifted));
    }

    private static LiftedAlignment rebuildLiftedWithCigar(final LiftedAlignment original, final String newCigar)
    {
        final LiftedAlignment rebuilt = new LiftedAlignment(
                original.Source, original.OrigContig, original.OrigPos, original.OrigCigar,
                original.LiftedChrom, original.LiftedPos, newCigar,
                original.AlignmentScore, original.NumMismatches,
                original.TransName, original.GeneId, original.GeneName,
                original.SoftClipAtBoundary, original.ForwardStrand);
        rebuilt.Dropped = original.Dropped;
        rebuilt.IsPrimaryChoice = original.IsPrimaryChoice;
        return rebuilt;
    }

    private LiftBackResult supplementaryResult(final SAMRecord record)
    {
        LiftedAlignment lifted = liftAlignment(
                LiftedAlignment.AlignmentSource.SELF,
                record.getReferenceName(), record.getAlignmentStart(), record.getCigarString(),
                getInt(record, AS_TAG), getInt(record, NM_TAG),
                !record.getReadNegativeStrandFlag());

        if(lifted == null)
            return unliftableResult(LiftBackResult.RecordRole.SUPPLEMENTARY, 0, "supp_translate_failed");

        if(lifted.fromTxContig())
        {
            final Cigar trimmed = ContigTranslator.trimTrailingMicroAnchor(
                    TextCigarCodec.decode(lifted.LiftedCigar), ANNOTATED_JUNCTION_MIN_ANCHOR_BP);
            if(!trimmed.toString().equals(lifted.LiftedCigar))
                lifted = rebuildLiftedWithCigar(lifted, trimmed.toString());
        }

        lifted.IsPrimaryChoice = true;

        // supplementaries don't carry XA, so the collapsed-locus check doesn't apply. Rescue MAPQ only on
        // a tx-contig supp where the input dropped to 0 due to the multi-alt-contig tie artefact.
        final int inputMapq = record.getMappingQuality();
        final int suppMapq = (lifted.fromTxContig() && inputMapq == 0) ? RESCUED_MAPQ : inputMapq;

        final int numRefAlts = lifted.fromTxContig() ? 0 : 1;
        final int numTxAlts = lifted.fromTxContig() ? 1 : 0;

        return new LiftBackResult(
                LiftBackCategory.SUPPLEMENTARY, LiftBackResult.Composition.fromAlignments(List.of(lifted)),
                LiftBackResult.RecordRole.SUPPLEMENTARY,
                lifted.LiftedChrom, lifted.LiftedPos, lifted.LiftedCigar,
                !lifted.ForwardStrand,
                lifted.cigarHasN(), inputMapq, suppMapq,
                0, numRefAlts, numTxAlts,
                1, 1,
                false, false, false, false,
                lifted.GeneId != null ? lifted.GeneId : "", "",
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

        return new LiftedAlignment(
                source, contig, pos, cigarStr,
                translated.chromosome(), translated.genomicStart(), translated.genomicCigar().toString(),
                as, nm,
                entry.transName(), entry.geneId(), entry.geneName(),
                softClipAtBoundary, forwardStrand);
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
                List.of());
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
                List.of());
    }
}
