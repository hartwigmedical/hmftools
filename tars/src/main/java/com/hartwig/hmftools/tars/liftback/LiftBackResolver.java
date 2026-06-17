package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.common.SpliceCommon.ALT_CONTIG_SUFFIX;
import static com.hartwig.hmftools.tars.common.TarsConfig.TARS_LOGGER;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.tars.common.BwaMemScore;
import com.hartwig.hmftools.tars.common.ContigEntry;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;

// Lifts SAMRecords aligned against ref + Tx contigs to genomic coords. Every input record produces exactly one
// result; UNMAPPED / LIFT_FAILED results carry placeholder fields so the emit step can flag records unmapped.
public class LiftBackResolver
{
    private static final String XA_TAG = "XA";
    private static final String AS_TAG = "AS";
    private static final String XS_TAG = "XS";
    private static final String NM_TAG = "NM";

    // bwa's unique-placement MAPQ, also the value liftback assigns when it resolves a placement.
    private static final int RESCUE_MAPQ = 60;

    // Set to 1: every tx-contig junction is annotated by construction, so even a 1bp anchor is a real exon
    // base. A higher floor manufactured false junction diffs by clipping legitimate short anchors.
    private static final int ANNOTATED_JUNCTION_MIN_ANCHOR_BP = 1;

    // Higher than bare floor because an adjacent softclip indicates bwa over-ran the exon boundary -
    // the tiny anchor is an over-extension artefact, not the read's real terminus.
    private static final int ANNOTATED_JUNCTION_MIN_SOFTCLIP_ANCHOR_BP = 3;

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
            mSegmentsByAltContig.computeIfAbsent(entry.contigName(), k -> new ArrayList<>()).add(entry);
        for(final List<ContigEntry> segments : mSegmentsByAltContig.values())
            segments.sort(Comparator.comparingInt(ContigEntry::altStart));

        mExonIndex = exonIndex;
    }

    public Set<String> contigNames()
    {
        return mSegmentsByAltContig.keySet();
    }

    // Returns the segment owning altPos, or null if it falls in an inter-transcript spacer or the contig is unknown.
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
            // altPos is before the first segment; return it so ContigTranslator's leading-overhang clamp can salvage it.
            return segments.isEmpty() ? null : segments.get(0);
        }

        final ContigEntry segment = segments.get(candidate);
        if(altPos <= segment.altEnd())
            return segment;

        // altPos in the inter-segment spacer; choose the nearer neighbour and let ContigTranslator clamp the overhang.
        if(candidate + 1 < segments.size())
        {
            final ContigEntry next = segments.get(candidate + 1);
            final int leadingOverhang = next.altStart() - altPos;
            final int trailingOverhang = altPos - segment.altEnd();
            return leadingOverhang <= trailingOverhang ? next : segment;
        }

        return segment;
    }

    // Lift-only entry point for callers that don't need the full LiftBackResult (e.g. SA tag rewriting).
    public LiftedCoords liftCoords(final String contig, final int pos, final String cigarStr)
    {
        final LiftedAlignment lifted = liftAlignment(
                LiftedAlignment.AlignmentSource.SELF, contig, pos, cigarStr, 0, 0, true);
        if(lifted == null)
            return null;
        return new LiftedCoords(lifted.LiftedChrom, lifted.LiftedPos, lifted.LiftedCigar);
    }

    private LiftedAlignment liftSelf(final SAMRecord record)
    {
        return liftAlignment(
                LiftedAlignment.AlignmentSource.SELF,
                record.getReferenceName(), record.getAlignmentStart(), record.getCigarString(),
                getInt(record, AS_TAG), getInt(record, NM_TAG), !record.getReadNegativeStrandFlag());
    }

    public LiftBackResult resolve(final SAMRecord record)
    {
        if(record.getReadUnmappedFlag())
            return unmappedResult(record);

        if(record.getSupplementaryAlignmentFlag())
            return supplementaryResult(record);

        return resolvePrimary(record);
    }

    private LiftBackResult resolvePrimary(final SAMRecord record)
    {
        final List<LiftedAlignment> xaAlts = parseAndLiftXa(record);
        return resolvePrimaryWithAlts(record, xaAlts, xaAlts.size());
    }

    private LiftBackResult resolvePrimaryWithAlts(
            final SAMRecord record, final List<LiftedAlignment> alts, final int numXaAltsForReport)
    {
        final LiftedAlignment self = liftSelf(record);

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
                .toList();

        final int numLoci = countDistinctLoci(keptAlignments);
        final int cigarsAtPrimaryLocus = countDistinctCigarsAtLocus(keptAlignments, effectivePrimary);
        final String geneIds = joinGeneIds(keptAlignments);

        final int inputMapq = record.getMappingQuality();
        final boolean swapped = effectivePrimary != self;
        final boolean hiddenTie = inputMapq == 0 && hasHiddenTie(record);
        final boolean inAnnotatedExon = mExonIndex != null
                && mExonIndex.contains(effectivePrimary.LiftedChrom, effectivePrimary.LiftedPos);
        final boolean hasTxMatch = features.NumTxAlts > 0;
        final int updatedMapq = decidePrimaryMapq(
                inputMapq, numLoci, swapped, hiddenTie, effectivePrimary.fromTxContig(), inAnnotatedExon, hasTxMatch);

        final LiftedAlignment primaryCoords = effectivePrimary;

        if(swapped)
            TARS_LOGGER.trace("discriminator {} {}: primary -> {}:{} {} ({})",
                    features.Category, record.getReadName(), primaryCoords.LiftedChrom, primaryCoords.LiftedPos,
                    primaryCoords.LiftedCigar, outcome.note());

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

    // MAPQ=0 on a tx-contig supplementary is the multi-alt-contig tie artefact; rescue to RESCUE_MAPQ.
    private LiftBackResult supplementaryResult(final SAMRecord record)
    {
        final LiftedAlignment lifted = liftSelf(record);

        if(lifted == null)
            return unliftableResult(LiftBackResult.RecordRole.SUPPLEMENTARY, 0, "supp_translate_failed");

        lifted.IsPrimaryChoice = true;

        final int inputMapq = record.getMappingQuality();
        final int outputMapq = (lifted.fromTxContig() && inputMapq == 0) ? RESCUE_MAPQ : inputMapq;
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

    // Self is excluded from the dedup key set so a Tx XA alt lifting to the same coords as a ref self is preserved (drives BOTH_AGREE).
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

            // bwa XA pos is signed (sign encodes strand). Tolerate malformed entries.
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

        // Drop sub-threshold terminal anchors that bwa fabricated by walking past an exon boundary.
        final ContigTranslator.MicroAnchorResult trimmed = ContigTranslator.trimMicroAnchors(
                translated.genomicCigar(), ANNOTATED_JUNCTION_MIN_ANCHOR_BP, ANNOTATED_JUNCTION_MIN_SOFTCLIP_ANCHOR_BP);

        return new LiftedAlignment(
                source, contig, pos, cigarStr,
                translated.chromosome(), translated.genomicStart() + trimmed.StartShift,
                trimmed.AdjustedCigar.toString(),
                as, nm,
                entry.transName(), entry.geneId(), entry.geneName(),
                softClipAtBoundary, forwardStrand, entry.strand());
    }

    // Counts distinct loci among only the best-scoring alignments (sub-optimal alts excluded from the multimap count).
    // Sub-optimal XA alts aren't real placement competitors -- counting them inflated numLoci and blocked MAPQ rescue.
    private static int countDistinctLoci(final List<LiftedAlignment> alignments)
    {
        if(alignments.isEmpty())
            return 0;

        int bestScore = Integer.MIN_VALUE;
        for(final LiftedAlignment la : alignments)
            bestScore = Math.max(bestScore, reconstructedScore(la));

        final Set<String> loci = new HashSet<>();
        for(final LiftedAlignment la : alignments)
        {
            if(reconstructedScore(la) == bestScore)
                loci.add(locusKey(la));
        }
        return loci.size();
    }

    // Reconstruct bwa-mem2 AS from CIGAR + NM (NM = mismatches + indel bases, so mismatches = NM - indelBases).
    static int reconstructedScore(final LiftedAlignment alignment)
    {
        final Cigar cigar = TextCigarCodec.decode(alignment.OrigCigar);
        int matched = 0;
        int indelBases = 0;
        int gapOps = 0;
        for(final CigarElement element : cigar.getCigarElements())
        {
            final CigarOperator op = element.getOperator();
            if(op == CigarOperator.M || op == CigarOperator.EQ || op == CigarOperator.X)
                matched += element.getLength();
            else if(op == CigarOperator.I || op == CigarOperator.D)
            {
                indelBases += element.getLength();
                ++gapOps;
            }
        }

        // XA alts carry no AS tag, so score is reconstructed from the CIGAR to filter out sub-optimal competitors.
        final int mismatches = Math.max(0, alignment.NumMismatches - indelBases);
        return (matched - mismatches) * BwaMemScore.MATCH + mismatches * BwaMemScore.MISMATCH
                + gapOps * BwaMemScore.GAP_OPEN + indelBases * BwaMemScore.GAP_EXTEND;
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

    // Rescue to RESCUE_MAPQ if: swapped by discriminator; or MAPQ=0 + single locus + no unresolved hidden tie.
    // Hidden tie (XS==AS, ref-only primary, not in annotated exon) blocks rescue - the unseen alt may be real.
    // Gated on hasTxMatch: without a tx alignment, MAPQ-0 is not a tx-artefact and is left alone.
    static int decidePrimaryMapq(
            final int inputMapq, final int numLoci, final boolean swapped, final boolean hiddenTie,
            final boolean primaryFromTxContig, final boolean primaryInAnnotatedExon, final boolean hasTxMatch)
    {
        if(!hasTxMatch)
            return inputMapq;
        if(swapped)
            return RESCUE_MAPQ;
        final boolean unresolvedHiddenTie = hiddenTie && !primaryFromTxContig && !primaryInAnnotatedExon;
        if(numLoci == 1 && inputMapq == 0 && !unresolvedHiddenTie)
            return RESCUE_MAPQ;
        return inputMapq;
    }

    // When XS == AS, an equally-scoring alt wasn't emitted by bwa. Flag as a hidden tie to skip MAPQ rescue.
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
