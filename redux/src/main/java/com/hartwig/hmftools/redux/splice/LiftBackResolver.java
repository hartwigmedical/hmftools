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
    private static final String NM_TAG = "NM";

    private static final int RESCUED_MAPQ = 60;

    // minimum flanking M length for an N operator to be treated as a real splice junction in
    // discriminator decisions. Lift-back of a tx alignment whose last bases bleed into the next exon
    // can produce N with anchors of 1-3 bp — those are not evidence of splicing and must not swap the
    // primary off a clean ref full-match. 8 bp gives ~1-in-65k random-match probability.
    static final int MIN_REAL_N_ANCHOR = 8;

    // per-alt-contig list of segments sorted by altStart so a record's alt-contig position can be
    // bin-searched back to the owning transcript. Non-alt-contig alignments (ref) fall through unindexed.
    private final Map<String, List<ContigEntry>> mSegmentsByAltContig;

    public LiftBackResolver(final List<ContigEntry> entries)
    {
        mSegmentsByAltContig = new HashMap<>();
        for(ContigEntry entry : entries)
            mSegmentsByAltContig.computeIfAbsent(entry.contigName(), k -> new ArrayList<>()).add(entry);
        for(List<ContigEntry> segments : mSegmentsByAltContig.values())
            segments.sort(Comparator.comparingInt(ContigEntry::altStart));
    }

    public Set<String> contigNames()
    {
        return mSegmentsByAltContig.keySet();
    }

    // locates the transcript segment that owns altPos on the given alt contig, or null if the position
    // falls outside any transcript (e.g. inside the inter-transcript N spacer) or the contig isn't an alt contig.
    ContigEntry findSegment(final String altContig, final int altPos)
    {
        List<ContigEntry> segments = mSegmentsByAltContig.get(altContig);
        if(segments == null)
            return null;

        // binary search: find the rightmost segment with altStart <= altPos, then check altEnd
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
            // altPos sits before the first segment's altStart (in the upstream spacer). Hand back the first
            // segment so ContigTranslator's leading-overhang clamp can salvage it.
            return segments.isEmpty() ? null : segments.get(0);
        }

        ContigEntry segment = segments.get(candidate);
        if(altPos <= segment.altEnd())
            return segment;

        // altPos sits in the spacer between candidate and candidate+1. Choose whichever neighbour the read
        // overhangs less, and let ContigTranslator's clamp convert the overhang into soft-clip.
        if(candidate + 1 < segments.size())
        {
            ContigEntry next = segments.get(candidate + 1);
            int leadingOverhang = next.altStart() - altPos;
            int trailingOverhang = altPos - segment.altEnd();
            return leadingOverhang <= trailingOverhang ? next : segment;
        }

        return segment;
    }

    // translates a single (contig, pos, cigar) triple from transcript-contig coordinates to genomic.
    // pass-through (returns the input values) for non-alt contigs. Returns null when the contig is an alt
    // contig but the position cannot be translated (e.g. it falls into the inter-transcript N spacer).
    //
    // intended for callers that want lift-only translation without the full LiftBackResult machinery
    // (e.g. SA tag rewriting, where mapQ / NM pass through unchanged).
    public LiftedCoords liftCoords(final String contig, final int pos, final String cigarStr)
    {
        LiftCore core = lift(contig, pos, cigarStr);
        return core == null ? null : new LiftedCoords(core.LiftedChrom, core.LiftedPos, core.LiftedCigar);
    }

    public static class LiftedCoords
    {
        public final String Chromosome;
        public final int Position;
        public final String CigarString;

        public LiftedCoords(final String chromosome, final int position, final String cigarString)
        {
            Chromosome = chromosome;
            Position = position;
            CigarString = cigarString;
        }
    }

    // shared lift kernel for liftCoords + liftAlignment. Returns null when the contig is an alt contig but
    // the position cannot be translated; ref-genome alignments fall through with Entry / ParsedCigar = null.
    private static final class LiftCore
    {
        final ContigEntry Entry;            // null = ref pass-through
        final String LiftedChrom;
        final int LiftedPos;
        final String LiftedCigar;
        final Cigar ParsedCigar;            // null = ref pass-through (no exon-boundary check needed)

        LiftCore(final ContigEntry entry, final String liftedChrom, final int liftedPos, final String liftedCigar,
                final Cigar parsedCigar)
        {
            Entry = entry;
            LiftedChrom = liftedChrom;
            LiftedPos = liftedPos;
            LiftedCigar = liftedCigar;
            ParsedCigar = parsedCigar;
        }
    }

    private LiftCore lift(final String contig, final int pos, final String cigarStr)
    {
        if(!mSegmentsByAltContig.containsKey(contig))
        {
            // alt contig missing from the segments map (FASTA/sidecar mismatch) must not pass through
            // as if it were ref — the resulting "lifted" coords would leak _tx contig names into the BAM
            // (read RNAME, mate RNEXT, MC) and break downstream tools that key off genomic coords.
            if(contig.endsWith(ALT_CONTIG_SUFFIX))
                return null;

            // ref alignment — pass through unchanged
            return new LiftCore(null, contig, pos, cigarStr, null);
        }

        ContigEntry entry = findSegment(contig, pos);
        if(entry == null)
            return null;

        Cigar cigar = TextCigarCodec.decode(cigarStr);
        ContigTranslator.TranslationResult result = ContigTranslator.translate(entry, pos, cigar);
        if(result == null)
            return null;

        return new LiftCore(entry, result.chromosome(), result.genomicStart(), result.genomicCigar().toString(), cigar);
    }

    public LiftBackResult resolve(final SAMRecord record)
    {
        if(record.getReadUnmappedFlag())
            return unmappedResult(record);

        if(record.isSecondaryOrSupplementary())
            return supplementaryResult(record);

        return resolvePrimary(record);
    }

    private LiftBackResult resolvePrimary(final SAMRecord record)
    {
        LiftedAlignment selfAlignment = liftAlignment(
                LiftedAlignment.AlignmentSource.SELF,
                record.getReferenceName(), record.getAlignmentStart(), record.getCigarString(),
                getInt(record, AS_TAG), getInt(record, NM_TAG),
                !record.getReadNegativeStrandFlag());

        // bwa's chosen primary failed to translate — record kept but flagged unmapped on emit
        if(selfAlignment == null)
            return unliftableResult(LiftBackResult.RecordRole.PRIMARY, countXaEntries(record), "primary_translate_failed");

        selfAlignment.IsPrimaryChoice = true;

        List<LiftedAlignment> xaAlts = parseAndLiftXa(record);

        List<LiftedAlignment> allAlignments = new ArrayList<>(1 + xaAlts.size());
        allAlignments.add(selfAlignment);
        allAlignments.addAll(xaAlts);

        Features features = categorize(allAlignments);

        // confident discriminator outcomes may swap the BAM primary (when bwa picked a paralog over the
        // spliced transcript alignment) and mark losing-side alts with Dropped=true. Returns the effective
        // primary — usually selfAlignment, but a winning tx alt when a swap occurred.
        DiscriminatorOutcome outcome = applyDiscriminator(allAlignments, features.Category, selfAlignment);
        LiftedAlignment effectivePrimary = outcome.EffectivePrimary;

        // aggregate stats (loci, cigarsAtPrimary, geneIds) over the kept (non-Dropped) alignments,
        // matching what the BAM XA tag will surface.
        List<LiftedAlignment> keptAlignments = keptAlignments(allAlignments);
        Aggregates aggregates = aggregate(keptAlignments, effectivePrimary);

        // bwa drops MAPQ to 0 when a read ties across multiple alt transcript contigs of the same gene.
        // when every lifted alignment (self + XA alts) collapses to a single genomic locus, the tie is
        // an artefact of the transcript-contig representation, not a real multi-mapper — restore MAPQ
        // so downstream tools (Isofox, redux) don't discard the read as ambiguous.
        //
        // also rescue when we swapped the BAM primary to a different locus — bwa's MAPQ belonged to the
        // (rejected) original locus and shouldn't carry over.
        int bwaMapq = record.getMappingQuality();
        boolean swapped = effectivePrimary != selfAlignment;
        int updatedMapq;
        if(swapped)
            updatedMapq = RESCUED_MAPQ;
        else if(aggregates.NumLoci == 1 && bwaMapq == 0)
            updatedMapq = RESCUED_MAPQ;
        else
            updatedMapq = bwaMapq;

        return new LiftBackResult(
                features.Category, LiftBackResult.Composition.fromAlignments(keptAlignments),
                LiftBackResult.RecordRole.PRIMARY,
                effectivePrimary.LiftedChrom, effectivePrimary.LiftedPos, effectivePrimary.LiftedCigar,
                !effectivePrimary.ForwardStrand,
                effectivePrimary.cigarHasN(), bwaMapq, updatedMapq,
                xaAlts.size(), features.NumRefAlts, features.NumTxAlts,
                aggregates.NumLoci, aggregates.CigarsAtPrimaryLocus,
                features.TxHasNCigar, features.TxSoftClipAtBoundary,
                features.RefSoftClipped, features.RefFullMatch,
                aggregates.GeneIds, outcome.Note,
                allAlignments);
    }

    private static List<LiftedAlignment> keptAlignments(final List<LiftedAlignment> alignments)
    {
        return alignments.stream().filter(la -> !la.Dropped).collect(Collectors.toList());
    }

    private static class Aggregates
    {
        int NumLoci;
        int CigarsAtPrimaryLocus;
        String GeneIds;
    }

    private static Aggregates aggregate(final List<LiftedAlignment> alignments, final LiftedAlignment effectivePrimary)
    {
        String primaryLocus = locusKey(effectivePrimary);
        Set<String> loci = new HashSet<>();
        Set<String> cigarsAtPrimary = new HashSet<>();
        Set<String> geneIdSet = new HashSet<>();
        for(LiftedAlignment la : alignments)
        {
            String locus = locusKey(la);
            loci.add(locus);
            if(locus.equals(primaryLocus))
                cigarsAtPrimary.add(la.LiftedCigar);
            if(la.GeneId != null)
                geneIdSet.add(la.GeneId);
        }
        Aggregates aggregates = new Aggregates();
        aggregates.NumLoci = loci.size();
        aggregates.CigarsAtPrimaryLocus = cigarsAtPrimary.size();
        aggregates.GeneIds = geneIdSet.stream().sorted().collect(Collectors.joining("|"));
        return aggregates;
    }

    // bundles the effective-primary alignment with a diagnostic note. Effective primary equals selfAlignment
    // for most categories; for tx-favouring discriminators when bwa picked ref as primary, it is the winning
    // tx alt (BAM primary is swapped onto it, original ref-self demoted to an XA entry).
    private static class DiscriminatorOutcome
    {
        final LiftedAlignment EffectivePrimary;
        final String Note;

        DiscriminatorOutcome(final LiftedAlignment effectivePrimary, final String note)
        {
            EffectivePrimary = effectivePrimary;
            Note = note;
        }
    }

    // applies confident discriminator outcomes:
    //   - marks losing-side alts with Dropped=true (stay in LiftedAlignments for TSV-B diagnostics,
    //     excluded from the rebuilt BAM XA tag)
    //   - swaps the BAM primary onto the winning side when bwa's chosen primary disagreed with the
    //     discriminator (e.g. bwa picked a processed-pseudogene paralog over the spliced parent gene).
    private static DiscriminatorOutcome applyDiscriminator(
            final List<LiftedAlignment> alignments, final LiftBackCategory category, final LiftedAlignment self)
    {
        switch(category)
        {
            case BOTH_TX_JUNCTION_REF_SOFTCLIP:
            case BOTH_MULTI_TX_JUNCTION:
                return promoteTxOverRef(alignments, self);

            case BOTH_TX_JUNCTION_REF_MATCH:
                // ref full-match across the supposed intron with NM=0 is overwhelming evidence the read is
                // genuinely unspliced (pre-mRNA / retained-intron / DNA contamination). 100bp of intron
                // matching by chance is ~4^-100 — essentially zero. The tx N-CIGAR is the artifact, not the
                // ref full-match. Favour ref and drop the tx alt.
                return promoteRefOverTx(alignments, self);

            case BOTH_TX_SOFTCLIP_REF_MATCH:
                if(!self.fromTxContig())
                {
                    for(LiftedAlignment la : alignments)
                        if(la.fromTxContig())
                            la.Dropped = true;
                    return new DiscriminatorOutcome(self, "");
                }
                return new DiscriminatorOutcome(self, "self_was_tx_no_swap");

            default:
                return new DiscriminatorOutcome(self, "");
        }
    }

    // ref-favouring outcome (mirror of promoteTxOverRef). If self is ref, drop all tx alts. If self is
    // tx, swap: promote the winning ref alt to primary, demote self (tx) to XA, drop other tx alts.
    private static DiscriminatorOutcome promoteRefOverTx(
            final List<LiftedAlignment> alignments, final LiftedAlignment self)
    {
        if(!self.fromTxContig())
        {
            for(final LiftedAlignment la : alignments)
                if(la.fromTxContig())
                    la.Dropped = true;
            return new DiscriminatorOutcome(self, "");
        }

        // self is tx — promote the first ref alt
        LiftedAlignment winner = null;
        for(final LiftedAlignment la : alignments)
        {
            if(!la.fromTxContig())
            {
                winner = la;
                break;
            }
        }
        if(winner == null)
            return new DiscriminatorOutcome(self, "");

        self.IsPrimaryChoice = false;
        winner.IsPrimaryChoice = true;
        // keep self (tx) in XA as a transcript hit; drop OTHER tx alts
        for(final LiftedAlignment la : alignments)
        {
            if(la.fromTxContig() && la != self)
                la.Dropped = true;
        }
        return new DiscriminatorOutcome(winner, "swapped_tx_to_ref");
    }

    // tx-favouring outcome: keep the best tx alt as primary. If self is tx, drop all ref alts. If self is
    // ref, swap: promote the winning tx alt to primary, demote self to XA (preserves the original locus as
    // an informative paralog hit), and drop other ref alts.
    private static DiscriminatorOutcome promoteTxOverRef(
            final List<LiftedAlignment> alignments, final LiftedAlignment self)
    {
        if(self.fromTxContig())
        {
            for(LiftedAlignment la : alignments)
                if(!la.fromTxContig())
                    la.Dropped = true;
            return new DiscriminatorOutcome(self, "");
        }

        // self is ref — prefer the tx alt that actually carries the N junction; fall back to any tx alt
        LiftedAlignment winner = null;
        for(LiftedAlignment la : alignments)
        {
            if(la.fromTxContig() && la.cigarHasN())
            {
                winner = la;
                break;
            }
        }
        if(winner == null)
        {
            for(LiftedAlignment la : alignments)
            {
                if(la.fromTxContig())
                {
                    winner = la;
                    break;
                }
            }
        }
        if(winner == null)
            return new DiscriminatorOutcome(self, "");

        self.IsPrimaryChoice = false;
        winner.IsPrimaryChoice = true;
        // keep self in XA as a paralog reference; drop OTHER ref alts (paralog noise)
        for(LiftedAlignment la : alignments)
        {
            if(!la.fromTxContig() && la != self)
                la.Dropped = true;
        }
        return new DiscriminatorOutcome(winner, "swapped_ref_to_tx");
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

        lifted.IsPrimaryChoice = true;

        // supplementaries don't carry XA, so collapsed-locus check doesn't apply. Rescue only when bwa
        // dropped to 0 on a tx-contig supp (the multi-alt-contig tie artefact); leave non-zero MAPQ alone.
        int bwaMapq = record.getMappingQuality();
        int suppMapq = (lifted.fromTxContig() && bwaMapq == 0) ? RESCUED_MAPQ : bwaMapq;

        int numRefAlts = lifted.fromTxContig() ? 0 : 1;
        int numTxAlts = lifted.fromTxContig() ? 1 : 0;

        return new LiftBackResult(
                LiftBackCategory.SUPPLEMENTARY, LiftBackResult.Composition.fromAlignments(List.of(lifted)),
                LiftBackResult.RecordRole.SUPPLEMENTARY,
                lifted.LiftedChrom, lifted.LiftedPos, lifted.LiftedCigar,
                !lifted.ForwardStrand,
                lifted.cigarHasN(), bwaMapq, suppMapq,
                0, numRefAlts, numTxAlts,
                1, 1,
                false, false, false, false,
                lifted.GeneId != null ? lifted.GeneId : "", "",
                List.of(lifted));
    }

    // bundles the category decision with the per-alignment evidence the discriminator collected.
    // emitted into TSV-A so post-hoc analysis can answer "why did this read land in category X?"
    static class Features
    {
        LiftBackCategory Category;
        int NumRefAlts;
        int NumTxAlts;
        boolean TxHasNCigar;
        boolean TxSoftClipAtBoundary;
        boolean RefSoftClipped;
        boolean RefFullMatch;
        boolean RefHasNCigar;
    }

    private Features categorize(final List<LiftedAlignment> alignments)
    {
        Features features = new Features();

        if(alignments.isEmpty())
        {
            features.Category = LiftBackCategory.LIFT_FAILED;
            return features;
        }

        Set<String> loci = new HashSet<>();
        boolean anyHasN = false;
        Set<String> distinctCigars = new HashSet<>();

        for(LiftedAlignment alignment : alignments)
        {
            loci.add(locusKey(alignment));
            distinctCigars.add(alignment.LiftedCigar);
            if(alignment.cigarHasN())
                anyHasN = true;

            if(alignment.fromTxContig())
            {
                ++features.NumTxAlts;
                if(alignment.cigarHasRealNJunction(MIN_REAL_N_ANCHOR))
                    features.TxHasNCigar = true;
                if(alignment.SoftClipAtBoundary)
                    features.TxSoftClipAtBoundary = true;
            }
            else
            {
                ++features.NumRefAlts;
                if(alignment.cigarHasSoftClip())
                    features.RefSoftClipped = true;
                else
                    features.RefFullMatch = true;
                if(alignment.cigarHasRealNJunction(MIN_REAL_N_ANCHOR))
                    features.RefHasNCigar = true;
            }
        }

        boolean hasRef = features.NumRefAlts > 0;
        boolean hasTx = features.NumTxAlts > 0;

        if(loci.size() >= 2)
        {
            if(hasRef && hasTx)
            {
                // tx found an annotated junction (N-CIGAR) at one locus; ref alts hit other loci
                // intronlessly — almost always processed pseudogenes / paralogs. Favour tx and swap
                // the BAM primary onto it.
                if(features.TxHasNCigar && !features.RefHasNCigar)
                    features.Category = LiftBackCategory.BOTH_MULTI_TX_JUNCTION;
                else
                    features.Category = LiftBackCategory.BOTH_MULTI;
            }
            else if(hasTx)
                features.Category = LiftBackCategory.TX_MULTI;
            else
                features.Category = LiftBackCategory.REF_MULTI;
            return features;
        }

        // numLoci == 1
        if(hasRef && !hasTx)
        {
            features.Category = LiftBackCategory.REF_SINGLE;
            return features;
        }
        if(hasTx && !hasRef)
        {
            features.Category = LiftBackCategory.TX_SINGLE;
            return features;
        }

        // single locus, both ref and tx contributed — run discriminator
        if(distinctCigars.size() == 1 && !anyHasN)
            features.Category = LiftBackCategory.BOTH_AGREE;
        else if(features.TxHasNCigar && features.RefSoftClipped)
            features.Category = LiftBackCategory.BOTH_TX_JUNCTION_REF_SOFTCLIP;
        else if(features.TxHasNCigar && features.RefFullMatch && !features.RefSoftClipped)
            // tx found a real junction; ref ignored it and stretched through with a full match.
            // tx is the more faithful alignment — drop the ref alt downstream.
            features.Category = LiftBackCategory.BOTH_TX_JUNCTION_REF_MATCH;
        else if(features.TxSoftClipAtBoundary && features.RefFullMatch)
            features.Category = LiftBackCategory.BOTH_TX_SOFTCLIP_REF_MATCH;
        else
            features.Category = LiftBackCategory.BOTH_AMBIGUOUS;

        return features;
    }

    // parses the XA tag, lifts each alt, then dedupes (Stage 2):
    //   - drops XA entries that fail to lift or fail to parse
    //   - dedupes among XA alts by lifted (chrom, pos, CIGAR), keeping first occurrence.
    // Self is intentionally NOT included in the dedup key set so a Tx XA alt that lifts to the same
    // (chrom, pos, CIGAR) as a ref self-alignment is preserved (drives BOTH_AGREE).
    private List<LiftedAlignment> parseAndLiftXa(final SAMRecord record)
    {
        List<LiftedAlignment> alts = new ArrayList<>();
        String xa = record.getStringAttribute(XA_TAG);
        if(xa == null || xa.isEmpty())
            return alts;

        Set<String> seenKeys = new HashSet<>();

        for(String entry : xa.split(";"))
        {
            if(entry.isEmpty())
                continue;
            String[] parts = entry.split(",");
            if(parts.length < 4)
                continue;

            String contig = parts[0];
            String cigar = parts[2];

            // bwa XA pos is signed (sign = strand); be tolerant of malformed XA — skip entries we can't parse
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
                alts.add(lifted);
        }

        return alts;
    }

    private LiftedAlignment liftAlignment(
            final LiftedAlignment.AlignmentSource source, final String contig, final int pos, final String cigarStr,
            final int as, final int nm, final boolean forwardStrand)
    {
        LiftCore core = lift(contig, pos, cigarStr);
        if(core == null)
            return null;

        if(core.Entry == null)
        {
            // ref alignment — pass through unchanged; no transcript metadata, no boundary check
            return new LiftedAlignment(
                    source, contig, pos, cigarStr,
                    core.LiftedChrom, core.LiftedPos, core.LiftedCigar,
                    as, nm,
                    null, null, null,
                    false, forwardStrand);
        }

        boolean softClipAtBoundary = ContigTranslator.hasSoftClipAtExonBoundary(core.Entry, pos, core.ParsedCigar);

        return new LiftedAlignment(
                source, contig, pos, cigarStr,
                core.LiftedChrom, core.LiftedPos, core.LiftedCigar,
                as, nm,
                core.Entry.transName(), core.Entry.geneId(), core.Entry.geneName(),
                softClipAtBoundary, forwardStrand);
    }

    private static String liftedKey(final LiftedAlignment la)
    {
        return la.LiftedChrom + ":" + la.LiftedPos + ":" + la.LiftedCigar;
    }

    private static String locusKey(final LiftedAlignment la)
    {
        return la.LiftedChrom + ":" + la.LiftedPos;
    }

    private static int countXaEntries(final SAMRecord record)
    {
        String xa = record.getStringAttribute(XA_TAG);
        if(xa == null || xa.isEmpty())
            return 0;

        int count = 0;
        for(String entry : xa.split(";"))
        {
            if(!entry.isEmpty())
                ++count;
        }
        return count;
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
                LiftBackCategory.UNMAPPED, LiftBackResult.Composition.NONE,
                role,
                "*", 0, "*",
                false,
                false, 0, 0,
                0, 0, 0,
                0, 0,
                false, false, false, false,
                "", "",
                List.<LiftedAlignment>of());
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
                List.<LiftedAlignment>of());
    }
}
