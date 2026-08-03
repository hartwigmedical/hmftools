package com.hartwig.hmftools.tars.liftback.supplementary;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static com.hartwig.hmftools.common.utils.Arrays.reverseArray;

import static htsjdk.samtools.util.SequenceUtil.basesEqual;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.tars.common.BwaScoring;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

// Merges a primary's terminal softclip with an annotated-intron-spanning supplementary into a
// single spliced primary (N op). Chains up to SupplementaryConfig.MaxSuppMerges merges per read; README Step 3.
//
// A supplementary merges into the primary only when all of these hold. RejectReason names the first that fails; several
// reasons are the two-sided variants of one rule, so there are fewer conditions than reasons.
//   1. the primary has a terminal softclip to extend across      NO_TERMINAL_SOFTCLIP
//   2. exactly one supplementary reaches that side               MULTIPLE_SUPPS_IN_REACH
//   3. the supp is on the same chromosome and the same strand     DIFFERENT_CHROMOSOME, OPPOSITE_STRAND
//   4. the cigars are complementary and simple - no hard clip,
//      no indel at the boundary, full read length, one clear side  NO_MATCHING_SUPP, COMPLEX_CIGAR_SHAPE
//   5. primary M and supp M abut on the read, overlapping by at
//      most MaxSuppReadOverlap and not at all on the reference     READ_COVERAGE_GAP, READ_COVERAGE_OVERLAP
//   6. the implied intron is within [MinIntronLength, MaxIntronLength]  INTRON_TOO_SHORT, INTRON_TOO_LONG
//   7. both anchors are at least MinAnchorOverhang                 SHORT_ANCHOR
//   8. a junction position scores above Tier.NONE, or is annotated
//      when AnnotatedOnly is set                                   NOVEL_JUNCTION
// With no supplementary at all, tryRefVerifyOnly takes over and reports the REF_VERIFY_* reasons instead.
// TODO revisit this later to simplify
public class SupplementaryResolver
{
    // Plain-value input, testable without SAMRecord. readBases seeds the tie-break between equally scoring
    // junction positions; mateHintIntrons biases that choice toward the partner mate's resolved junctions.
    // NOTE: generated equals/hashCode compare readBases by array identity, so do not use as a Set/Map key.
    public record Candidate(
            String chromosome, boolean forwardStrand, int readLength, int primaryStart, String primaryCigar,
            List<Supplementary> supplementaries, byte[] readBases, List<ChrBaseRegion> mateHintIntrons)
    {
        public Candidate
        {
            if(mateHintIntrons == null)
            {
                mateHintIntrons = Collections.emptyList();
            }
        }

        public Candidate(
                final String chromosome, final boolean forwardStrand, final int readLength,
                final int primaryStart, final String primaryCigar, final List<Supplementary> supplementaries)
        {
            this(chromosome, forwardStrand, readLength, primaryStart, primaryCigar, supplementaries, null, Collections.emptyList());
        }
    }

    // One supplementary under consideration as a merge partner. index is its position in the caller's supplementary
    // list so the resolver can report which supps to drop without exposing SAMRecord.
    public record Supplementary(
            int index, String chromosome, boolean forwardStrand, int start, String cigar, int mapQuality)
    {
    }

    // On success, mergedCigar/mergedStart describe the new primary and droppedSupplementaryIndices lists the absorbed
    // supps. On failure, rejectReason carries the gate hit, and adjustedCigar/adjustedStart are set when the boundary
    // was snapped back to a splice motif without merging - null/-1 when the boundary was left alone.
    public record Result(
            boolean merged, String mergedCigar, int mergedStart, List<Integer> droppedSupplementaryIndices,
            List<ChrBaseRegion> introducedIntrons, int chainDepth, RejectReason rejectReason)
    {
        private static final Result NO_MERGE_NO_OP = new Result(
                false, null, -1, Collections.emptyList(), Collections.emptyList(), 0, RejectReason.NO_TERMINAL_SOFTCLIP);

        public static Result noMerge(final RejectReason reason)
        {
            if(reason == RejectReason.NO_TERMINAL_SOFTCLIP)
            {
                return NO_MERGE_NO_OP;
            }
            return new Result(false, null, -1, Collections.emptyList(), Collections.emptyList(), 0, reason);
        }
    }

    // Why a candidate was not merged, reported on the result.
    public enum RejectReason
    {
        NO_TERMINAL_SOFTCLIP,         // primary cigar has no leading or trailing S to extend across
        NO_MATCHING_SUPP,             // no supplementary had a complementary cigar shape
        DIFFERENT_CHROMOSOME,
        OPPOSITE_STRAND,
        READ_COVERAGE_OVERLAP,        // primary's M + supp's M overlap on the read
        READ_COVERAGE_GAP,            // primary's M + supp's M leave a gap on the read
        INTRON_TOO_SHORT,
        INTRON_TOO_LONG,
        SHORT_ANCHOR,                 // primary or supp matched portion < MinAnchorOverhang
        NOVEL_JUNCTION,               // candidate intron not in annotated set (and AnnotatedOnly is true)
        COMPLEX_CIGAR_SHAPE,          // hard clip, indel adjacent to softclip boundary, etc.
        MULTIPLE_SUPPS_IN_REACH,      // more than one supp within merge reach; refuse to guess which splice
        // ref-verify path: primary has a terminal softclip but no matching supp
        REF_VERIFY_NO_CANDIDATE_EXON, // no annotated junction's adjacent exon lines up with the softclip
        REF_VERIFY_MISMATCH_TOO_HIGH, // softclipped read bases don't match the candidate exon's ref bases
        REF_VERIFY_AMBIGUOUS          // multiple annotated downstream/upstream exons match; refuse to guess
    }

    // Splice-motif strength of a candidate junction, weakest to strongest. Declaration order is meaningful: candidates
    // are ranked with compareTo, so a stronger motif outranks a weaker one.
    public enum Tier
    {
        NONE,            // neither donor nor acceptor matches a known splice motif
        SEMI_CANONICAL,  // GC-AG / AT-AC (and reverse-complement equivalents)
        CANONICAL,       // GT-AG (~99% of splice sites)
        ANNOTATED        // matches an annotated junction from the sidecar
    }

    // no junction position scored above Tier.NONE; a real position is a read offset, so never negative
    private static final int NO_JUNCTION_POSITION = -1;

    private final AnnotatedJunctionIndex mAnnotatedIndex;
    private final RefGenomeInterface mRefGenome;
    private final SupplementaryConfig mConfig;

    public SupplementaryResolver(final Set<ChrBaseRegion> annotatedJunctions, final SupplementaryConfig config)
    {
        this(
                new AnnotatedJunctionIndex(annotatedJunctions != null ? annotatedJunctions : new HashSet<>()),
                null, config);
    }

    public SupplementaryResolver(
            final AnnotatedJunctionIndex annotatedIndex, final RefGenomeInterface refGenome,
            final SupplementaryConfig config)
    {
        mAnnotatedIndex = annotatedIndex != null ? annotatedIndex : new AnnotatedJunctionIndex(new HashSet<>());
        mRefGenome = refGenome;
        mConfig = config;
    }

    private byte[] refBases(final String chromosome, final int posStart, final int posEnd)
    {
        return BwaScoring.refWindow(mRefGenome, chromosome, posStart, posEnd);
    }

    public Result resolve(final Candidate candidate)
    {
        if(candidate.supplementaries().isEmpty())
        {
            return tryRefVerifyOnly(candidate);
        }

        List<CigarElement> primaryCigar = CigarUtils.cigarElementsFromStr(candidate.primaryCigar());
        if(CigarUtils.hasHardClip(primaryCigar))
        {
            return Result.noMerge(RejectReason.COMPLEX_CIGAR_SHAPE);
        }

        int primaryStart = candidate.primaryStart();
        List<Supplementary> remaining = new ArrayList<>(candidate.supplementaries());
        List<Integer> dropped = new ArrayList<>();
        List<ChrBaseRegion> introns = new ArrayList<>();
        RejectReason lastReject = null;
        int chainDepth = 0;

        while(chainDepth < mConfig.MaxSuppMerges && !remaining.isEmpty())
        {
            boolean primaryHasLeadingS = !primaryCigar.isEmpty()
                    && primaryCigar.get(0).getOperator() == CigarOperator.S;
            boolean primaryHasTrailingS = !primaryCigar.isEmpty()
                    && primaryCigar.get(primaryCigar.size() - 1).getOperator() == CigarOperator.S;
            if(!primaryHasLeadingS && !primaryHasTrailingS)
            {
                if(chainDepth == 0)
                {
                    return Result.noMerge(RejectReason.NO_TERMINAL_SOFTCLIP);
                }
                break;
            }

            MergeOutcome merge = pickBestSupplementary(
                    candidate, primaryStart, primaryCigar, remaining);

            if(!merge.isSuccess())
            {
                if(chainDepth == 0)
                {
                    return Result.noMerge(merge.Reject);
                }
                lastReject = merge.Reject;
                break;
            }

            primaryStart = merge.MergedStart;
            primaryCigar = merge.MergedCigar;
            dropped.add(merge.MergedSupp.index());
            introns.add(merge.IntroducedIntron);
            remaining.remove(merge.MergedSupp);
            ++chainDepth;
        }

        if(chainDepth == 0)
        {
            // Belt and braces: initial-iteration failures should have returned above.
            return Result.noMerge(lastReject != null ? lastReject : RejectReason.NO_MATCHING_SUPP);
        }

        return new Result(
                true, CigarUtils.cigarElementsToStr(primaryCigar), primaryStart,
                dropped, introns, chainDepth, null);
    }

    private Result tryRefVerifyOnly(final Candidate candidate)
    {
        if(mRefGenome == null || candidate.readBases() == null)
        {
            return Result.noMerge(RejectReason.NO_MATCHING_SUPP);
        }

        List<CigarElement> primaryCigar = CigarUtils.cigarElementsFromStr(candidate.primaryCigar());
        if(CigarUtils.hasHardClip(primaryCigar))
        {
            return Result.noMerge(RejectReason.COMPLEX_CIGAR_SHAPE);
        }
        if(primaryCigar.isEmpty())
        {
            return Result.noMerge(RejectReason.NO_TERMINAL_SOFTCLIP);
        }

        boolean trailingS = primaryCigar.get(primaryCigar.size() - 1).getOperator() == CigarOperator.S;
        boolean leadingS = primaryCigar.get(0).getOperator() == CigarOperator.S;
        if(!trailingS && !leadingS)
        {
            return Result.noMerge(RejectReason.NO_TERMINAL_SOFTCLIP);
        }

        // Both-end clips: try longer clip first (more likely the splice tail). First win is taken.
        boolean[] sides;
        if(trailingS && leadingS)
        {
            boolean trailingFirst =
                    CigarUtils.rightSoftClipLength(primaryCigar) >= CigarUtils.leftSoftClipLength(primaryCigar);
            sides = new boolean[] { trailingFirst, !trailingFirst };
        }
        else
        {
            sides = new boolean[] { trailingS };
        }

        boolean anyCandidate = false;
        Result lastFailure = Result.noMerge(RejectReason.REF_VERIFY_MISMATCH_TOO_HIGH);
        for(boolean rightExtend : sides)
        {
            Result sideResult = attemptRefVerifySide(candidate, primaryCigar, rightExtend);
            if(sideResult.merged())
            {
                return sideResult;
            }
            if(sideResult.rejectReason() != RejectReason.REF_VERIFY_NO_CANDIDATE_EXON)
            {
                anyCandidate = true;
                lastFailure = sideResult;
            }
        }

        if(!anyCandidate)
        {
            return Result.noMerge(RejectReason.REF_VERIFY_NO_CANDIDATE_EXON);
        }
        return lastFailure;
    }

    // Ref-verify one terminal softclip against annotated donors/acceptors with boundary snap.
    // Does NOT count the no-candidate reject; tryRefVerifyOnly aggregates across both ends.
    private Result attemptRefVerifySide(
            final Candidate candidate, final List<CigarElement> primaryCigar, final boolean rightExtend)
    {
        if(MergeCigarOps.opAdjacentToSoftClip(primaryCigar, !rightExtend))
        {
            return Result.noMerge(RejectReason.COMPLEX_CIGAR_SHAPE);
        }

        int softclipLen = rightExtend
                ? CigarUtils.rightSoftClipLength(primaryCigar)
                : CigarUtils.leftSoftClipLength(primaryCigar);
        int primaryAnchor = rightExtend
                ? MergeCigarOps.matchedRun(primaryCigar, true)
                : MergeCigarOps.matchedRun(primaryCigar, false);

        // BWA often over-extends a few bases past the true exon boundary. Snap back up to
        // MaxRefVerifyBoundaryShift (smallest shift first), rolling over-extension into the softclip.
        int primaryRefEnd = candidate.primaryStart() + CigarUtils.cigarAlignedLength(primaryCigar) - 1;
        boolean anyCandidate = false;
        Result lastFailure = Result.noMerge(RejectReason.REF_VERIFY_MISMATCH_TOO_HIGH);

        for(int shift = 0; shift <= mConfig.MaxRefVerifyBoundaryShift; ++shift)
        {
            if(primaryAnchor - shift < 1)
                break;
            List<CigarElement> shiftedCigar = shift == 0
                    ? primaryCigar
                    : MergeCigarOps.shiftBoundaryIntoSoftclip(primaryCigar, shift, rightExtend);
            if(shiftedCigar == null)
                break;

            int boundary = rightExtend ? (primaryRefEnd + 1 - shift) : (candidate.primaryStart() - 1 + shift);
            List<ChrBaseRegion> candidates = rightExtend
                    ? mAnnotatedIndex.introByStart(candidate.chromosome(), boundary)
                    : mAnnotatedIndex.introByEnd(candidate.chromosome(), boundary);
            if(candidates.isEmpty())
                continue;

            anyCandidate = true;
            Result result = verifyAgainstCandidates(
                    candidate, shiftedCigar, candidates, softclipLen + shift, rightExtend);
            if(result.merged())
            {
                return result;
            }
            lastFailure = result;
        }

        if(!anyCandidate)
        {
            return Result.noMerge(RejectReason.REF_VERIFY_NO_CANDIDATE_EXON);
        }
        return lastFailure;
    }

    private Result verifyAgainstCandidates(
            final Candidate candidate, final List<CigarElement> primaryCigar,
            final List<ChrBaseRegion> candidates, final int softclipLen, final boolean rightExtend)
    {
        if(candidates.isEmpty())
        {
            return Result.noMerge(RejectReason.REF_VERIFY_NO_CANDIDATE_EXON);
        }

        byte[] readBases = candidate.readBases();
        byte[] softclipBases = new byte[softclipLen];
        if(rightExtend)
        {
            System.arraycopy(readBases, readBases.length - softclipLen, softclipBases, 0, softclipLen);
        }
        else
        {
            System.arraycopy(readBases, 0, softclipBases, 0, softclipLen);
        }

        // Match from the junction-proximal end (bwa-mem score walk). The outer softclip residual stays
        // soft-clipped since bwa often carries adapter/low-quality bases there (e.g. 19S132M ->
        // 6S15M..N..130M, not 21M..N..130M). Proximal end is the high index for a leading softclip.
        boolean proximalAtEnd = !rightExtend;
        ChrBaseRegion chosen = null;
        int chosenRun = 0;
        int chosenMismatches = Integer.MAX_VALUE;
        boolean ambiguous = false;

        for(ChrBaseRegion candidateIntron : candidates)
        {
            int intronLength = candidateIntron.end() - candidateIntron.start() + 1;
            if(intronLength < mConfig.MinIntronLength || intronLength > mConfig.MaxIntronLength)
                continue;

            int refStart;
            int refEnd;
            if(rightExtend)
            {
                refStart = candidateIntron.end() + 1;
                refEnd = candidateIntron.end() + softclipLen;
            }
            else
            {
                refStart = candidateIntron.start() - softclipLen;
                refEnd = candidateIntron.start() - 1;
            }

            byte[] refBases = refBases(candidate.chromosome(), refStart, refEnd);
            if(refBases == null || refBases.length != softclipLen)
                continue;

            int run = proximalScoringRun(softclipBases, refBases, rightExtend);
            if(run == 0)
                continue;
            int mismatches = mismatchesInRun(softclipBases, refBases, softclipLen, run, proximalAtEnd);

            // Longest run wins; tie-break fewest mismatches. Two introns tying both -> ambiguous -> reject.
            if(run > chosenRun || (run == chosenRun && mismatches < chosenMismatches))
            {
                chosen = candidateIntron;
                chosenRun = run;
                chosenMismatches = mismatches;
                ambiguous = false;
            }
            else if(run == chosenRun && mismatches == chosenMismatches && !candidateIntron.equals(chosen))
            {
                ambiguous = true;
            }
        }

        if(chosen == null)
        {
            return Result.noMerge(RejectReason.REF_VERIFY_MISMATCH_TOO_HIGH);
        }
        if(ambiguous)
        {
            return Result.noMerge(RejectReason.REF_VERIFY_AMBIGUOUS);
        }

        return buildRefVerifyMerge(candidate, primaryCigar, chosen, chosenRun, softclipLen, rightExtend);
    }

    private Result buildRefVerifyMerge(
            final Candidate candidate, final List<CigarElement> primaryCigar,
            final ChrBaseRegion chosen, final int matchedRun, final int softclipLen, final boolean rightExtend)
    {
        int intronLength = chosen.end() - chosen.start() + 1;
        int residualSoftclip = softclipLen - matchedRun;
        List<CigarElement> merged = new ArrayList<>(primaryCigar.size() + 3);
        if(rightExtend)
        {
            for(int i = 0; i < primaryCigar.size() - 1; ++i)
            {
                merged.add(primaryCigar.get(i));
            }
            merged.add(new CigarElement(intronLength, CigarOperator.N));
            merged.add(new CigarElement(matchedRun, CigarOperator.M));
            if(residualSoftclip > 0)
            {
                merged.add(new CigarElement(residualSoftclip, CigarOperator.S));
            }
        }
        else
        {
            if(residualSoftclip > 0)
            {
                merged.add(new CigarElement(residualSoftclip, CigarOperator.S));
            }
            merged.add(new CigarElement(matchedRun, CigarOperator.M));
            merged.add(new CigarElement(intronLength, CigarOperator.N));
            for(int i = 1; i < primaryCigar.size(); ++i)
            {
                merged.add(primaryCigar.get(i));
            }
        }

        int mergedStart = rightExtend ? candidate.primaryStart() : (chosen.start() - matchedRun);
        return new Result(
                true, CigarUtils.cigarElementsToStr(merged), mergedStart,
                Collections.emptyList(), Collections.singletonList(chosen),
                1, null);
    }

    // Longest run matching ref from the junction-proximal end, scored with the shared bwa-mem model
    // (match +1, mismatch -4) so ref-verify matches collapse/tail-extend. Leading clips
    // reverse both arrays so the walk runs from the high-index proximal end outward.
    private static int proximalScoringRun(final byte[] softclipBases, final byte[] refBases, final boolean rightExtend)
    {
        if(rightExtend)
        {
            return BwaScoring.maxScoringPrefix(softclipBases, refBases);
        }
        return BwaScoring.maxScoringPrefix(reverseArray(softclipBases), reverseArray(refBases));
    }

    private static int mismatchesInRun(
            final byte[] softclipBases, final byte[] refBases, final int len, final int run, final boolean proximalAtEnd)
    {
        int mismatches = 0;
        for(int k = 1; k <= run; ++k)
        {
            int index = proximalAtEnd ? (len - k) : (k - 1);
            if(!basesEqual(softclipBases[index], refBases[index]))
            {
                ++mismatches;
            }
        }
        return mismatches;
    }

    private Tier classifyJunctionTier(final ChrBaseRegion candidateIntron)
    {
        if(mAnnotatedIndex.contains(candidateIntron))
        {
            return Tier.ANNOTATED;
        }
        if(mRefGenome == null)
        {
            return Tier.NONE;
        }
        byte[] donor = refBases(
                candidateIntron.Chromosome, candidateIntron.start(), candidateIntron.start() + 1);
        byte[] acceptor = refBases(
                candidateIntron.Chromosome, candidateIntron.end() - 1, candidateIntron.end());
        return motifTier(donor, acceptor);
    }

    // Donor/acceptor 2-base flanks: GT-AG canonical (~99% of sites), GC-AG and AT-AC semi-canonical. The strand is
    // unknown at scan time, so each motif's reverse complement is accepted too.
    static Tier motifTier(final byte[] donorBases, final byte[] acceptorBases)
    {
        if(donorBases == null || donorBases.length != 2 || acceptorBases == null || acceptorBases.length != 2)
        {
            return Tier.NONE;
        }

        String donor = upperCase(donorBases);
        String acceptor = upperCase(acceptorBases);

        if(donor.equals("GT") && acceptor.equals("AG") || donor.equals("CT") && acceptor.equals("AC"))
        {
            return Tier.CANONICAL;
        }
        if(donor.equals("GC") && acceptor.equals("AG") || donor.equals("CT") && acceptor.equals("GC"))
        {
            return Tier.SEMI_CANONICAL;
        }
        if(donor.equals("AT") && acceptor.equals("AC") || donor.equals("GT") && acceptor.equals("AT"))
        {
            return Tier.SEMI_CANONICAL;
        }
        return Tier.NONE;
    }

    // the reference is soft-masked, so repeat-region bases arrive lower case
    private static String upperCase(final byte[] bases)
    {
        return new String(bases, StandardCharsets.US_ASCII).toUpperCase();
    }

    private MergeOutcome pickBestSupplementary(
            final Candidate candidate, final int primaryStart,
            final List<CigarElement> primaryCigar, final List<Supplementary> supps)
    {
        MergeOutcome chosen = null;
        RejectReason lastReject = null;
        boolean rightResolved = false;
        boolean leftResolved = false;

        for(Supplementary supp : supps)
        {
            MergeOutcome outcome = tryMerge(candidate, primaryStart, primaryCigar, supp);
            if(!outcome.isSuccess())
            {
                lastReject = outcome.Reject;
                continue;
            }

            // One supp per terminal softclip is the most that can be resolved, so a second reaching the same side
            // means there is no telling which splice is real - give up on the whole read rather than guess.
            if(outcome.RightExtend ? rightResolved : leftResolved)
            {
                return MergeOutcome.reject(RejectReason.MULTIPLE_SUPPS_IN_REACH);
            }

            rightResolved |= outcome.RightExtend;
            leftResolved |= !outcome.RightExtend;

            if(chosen == null || isBetterMerge(outcome, chosen))
            {
                chosen = outcome;
            }
        }

        return chosen != null
                ? chosen
                : MergeOutcome.reject(lastReject != null ? lastReject : RejectReason.NO_MATCHING_SUPP);
    }

    // Higher MAPQ wins, then the smaller intron. Only ever compares one right-extend against one left-extend, since a
    // second supp on either side is rejected above, and the chain loop picks up the loser next pass.
    private static boolean isBetterMerge(final MergeOutcome outcome, final MergeOutcome chosen)
    {
        if(outcome.MergedSupp.mapQuality() != chosen.MergedSupp.mapQuality())
        {
            return outcome.MergedSupp.mapQuality() > chosen.MergedSupp.mapQuality();
        }
        return outcome.IntroducedIntron.baseLength() < chosen.IntroducedIntron.baseLength();
    }

    private MergeOutcome tryMerge(
            final Candidate candidate, final int primaryStart,
            final List<CigarElement> primaryCigar, final Supplementary supp)
    {
        if(!candidate.chromosome().equals(supp.chromosome()))
        {
            return MergeOutcome.reject(RejectReason.DIFFERENT_CHROMOSOME);
        }

        if(candidate.forwardStrand() != supp.forwardStrand())
        {
            return MergeOutcome.reject(RejectReason.OPPOSITE_STRAND);
        }

        List<CigarElement> suppCigar = CigarUtils.cigarElementsFromStr(supp.cigar());
        if(CigarUtils.hasHardClip(suppCigar))
        {
            return MergeOutcome.reject(RejectReason.COMPLEX_CIGAR_SHAPE);
        }

        // ContigTranslator can expand a cross-exon M into M-N-M. If the post-N M coincidentally
        // overlaps the primary's span, clamp the supp to its primary-distal anchor before merge logic.
        Side primarySide = Side.of(primaryStart, primaryCigar);
        int suppStart = supp.start();
        MergeCigarOps.ClampedSupp clamped = MergeCigarOps.clampSuppToPrimaryBoundary(suppCigar, suppStart, primaryStart, primarySide.RefEnd);
        if(clamped != null)
        {
            suppCigar = clamped.Cigar;
            suppStart = clamped.Start;
        }

        Side suppSide = Side.of(suppStart, suppCigar);
        boolean rightExtend = primarySide.TrailingS > 0 && suppSide.LeadingS > 0;
        boolean leftExtend = primarySide.LeadingS > 0 && suppSide.TrailingS > 0;

        if(!rightExtend && !leftExtend)
        {
            return MergeOutcome.reject(RejectReason.NO_MATCHING_SUPP);
        }

        // Clipped at both ends, so the supp could extend either way. Whichever side of the primary it actually sits
        // on settles it; if it sits cleanly on neither, the shape is not one to guess at.
        if(rightExtend && leftExtend)
        {
            if(suppSide.Start > primarySide.RefEnd && suppSide.RefEnd >= primarySide.Start)
            {
                leftExtend = false;
            }
            else if(suppSide.RefEnd < primarySide.Start && suppSide.Start <= primarySide.RefEnd)
            {
                rightExtend = false;
            }
            else
            {
                return MergeOutcome.reject(RejectReason.COMPLEX_CIGAR_SHAPE);
            }
        }

        boolean primaryIsUpstream = rightExtend;
        Side up = primaryIsUpstream ? primarySide : suppSide;
        Side down = primaryIsUpstream ? suppSide : primarySide;

        return mergeJunction(candidate, up, down, primaryIsUpstream, supp);
    }

    // Conditions 4-6 of the merge policy for a resolved up/down anchor pair, in the order they are cheapest to
    // establish. Returns the first failure, or null when the pair may merge.
    private RejectReason anchorPairReject(
            final Candidate candidate, final Side up, final Side down, final int overlap, final int intronLength)
    {
        if(CigarUtils.cigarBaseLength(up.Cigar) != candidate.readLength()
                || CigarUtils.cigarBaseLength(down.Cigar) != candidate.readLength())
        {
            return RejectReason.COMPLEX_CIGAR_SHAPE;
        }

        if(MergeCigarOps.opAdjacentToSoftClip(up.Cigar, false) || MergeCigarOps.opAdjacentToSoftClip(down.Cigar, true))
        {
            return RejectReason.COMPLEX_CIGAR_SHAPE;
        }

        if(overlap < 0)
        {
            return RejectReason.READ_COVERAGE_GAP;
        }

        if(overlap > mConfig.MaxSuppReadOverlap || down.Start <= up.RefEnd)
        {
            return RejectReason.READ_COVERAGE_OVERLAP;
        }

        if(intronLength < mConfig.MinIntronLength)
        {
            return RejectReason.INTRON_TOO_SHORT;
        }

        if(intronLength > mConfig.MaxIntronLength)
        {
            return RejectReason.INTRON_TOO_LONG;
        }

        return null;
    }

    // Direction-agnostic merge: validates the anchor pair, scores candidate junction positions by tier, falls
    // back to the mate hint then the midpoint. supp is threaded as result payload only.
    private MergeOutcome mergeJunction(
            final Candidate candidate, final Side up, final Side down,
            final boolean primaryIsUpstream, final Supplementary supp)
    {
        int upMatchedRead = candidate.readLength() - up.TrailingS;
        int overlap = upMatchedRead - down.LeadingS;

        // Intron length is invariant under the junction position - computed once, checked in the gate below.
        int intronLength = (down.Start - 1 - up.RefEnd) + overlap;

        RejectReason gateReject = anchorPairReject(candidate, up, down, overlap, intronLength);
        if(gateReject != null)
        {
            return MergeOutcome.reject(gateReject);
        }

        // Junction position priority: annotated boundary (sidecar) > canonical then semi-canonical motif; ties within
        // the tier broken pseudo-randomly but deterministically (read-seeded) so ambiguous junctions distribute yet
        // stay reproducible. Falls back to the mate's intron, then the midpoint of the ambiguous range.
        int junctionReadPosition = scanJunctionPositions(candidate, up, down, upMatchedRead, primaryIsUpstream);

        if(junctionReadPosition == NO_JUNCTION_POSITION)
        {
            junctionReadPosition = junctionPositionFromMate(candidate, up, down, upMatchedRead, primaryIsUpstream);
        }

        if(junctionReadPosition == NO_JUNCTION_POSITION)
        {
            if(mConfig.AnnotatedOnly)
            {
                return MergeOutcome.reject(RejectReason.NOVEL_JUNCTION);
            }

            // midpoint of the ambiguous read range [down.LeadingS, upMatchedRead], rounded down. The scan and mate
            // paths check the anchors before accepting a position, so only this fallback has to.
            junctionReadPosition = (upMatchedRead + down.LeadingS) / 2;
            if(up.TrailingM < upMatchedRead - junctionReadPosition
                    || down.LeadingM < junctionReadPosition - down.LeadingS)
            {
                return MergeOutcome.reject(RejectReason.SHORT_ANCHOR);
            }
        }

        int upLoss = upMatchedRead - junctionReadPosition;
        int downLoss = junctionReadPosition - down.LeadingS;
        List<CigarElement> merged = MergeCigarOps.buildMergedCigar(up.Cigar, down.Cigar, upLoss, downLoss, intronLength);
        return MergeOutcome.success(
                up.Start, merged, intronAt(candidate, up, down, upMatchedRead, junctionReadPosition),
                supp, primaryIsUpstream);
    }

    // The intron implied by putting the junction at this read position: the upstream alignment gives up the bases
    // past it, the downstream alignment the bases before it.
    private static ChrBaseRegion intronAt(
            final Candidate candidate, final Side up, final Side down,
            final int upMatchedRead, final int readPosition)
    {
        int upLoss = upMatchedRead - readPosition;
        int downLoss = readPosition - down.LeadingS;
        return new ChrBaseRegion(
                candidate.chromosome(), up.RefEnd - upLoss + 1, down.Start + downLoss - 1);
    }

    // Scores every valid junction position by tier and returns one at the highest tier present. Ties at that tier are
    // broken pseudo-randomly but deterministically (read-seeded) so ambiguous junctions distribute yet stay
    // reproducible. NO_JUNCTION_POSITION when no position reaches a splice motif or annotated boundary.
    private int scanJunctionPositions(
            final Candidate candidate, final Side up, final Side down,
            final int upMatchedRead, final boolean primaryIsUpstream)
    {
        Tier bestTier = Tier.NONE;
        List<Integer> bestPositions = new ArrayList<>();

        int firstPosition = primaryIsUpstream ? upMatchedRead : down.LeadingS;
        int lastPosition = primaryIsUpstream ? down.LeadingS : upMatchedRead;
        int step = primaryIsUpstream ? -1 : 1;
        for(int readPosition = firstPosition;
                primaryIsUpstream ? readPosition >= lastPosition : readPosition <= lastPosition;
                readPosition += step)
        {
            if(up.TrailingM < upMatchedRead - readPosition || down.LeadingM < readPosition - down.LeadingS)
                continue;

            Tier tier = classifyJunctionTier(intronAt(candidate, up, down, upMatchedRead, readPosition));
            if(tier == Tier.NONE)
                continue;

            int cmp = tier.compareTo(bestTier);
            if(cmp > 0)
            {
                bestTier = tier;
                bestPositions.clear();
            }
            if(cmp >= 0)
            {
                bestPositions.add(readPosition);
            }
        }

        if(bestPositions.isEmpty())
        {
            return NO_JUNCTION_POSITION;
        }

        int index = bestPositions.size() == 1 ? 0 : Math.floorMod(tieBreakSeed(candidate, up), bestPositions.size());
        return bestPositions.get(index);
    }

    // Deterministic per-read seed so an ambiguous tie is reproducible run to run: hashes the read bases (when
    // available) with the upstream boundary so two reads over the same junction can pick differently.
    private static int tieBreakSeed(final Candidate candidate, final Side up)
    {
        int base = candidate.readBases() != null ? Arrays.hashCode(candidate.readBases()) : 0;
        return 31 * base + up.RefEnd;
    }

    // Use the mate's already-resolved intron as a hint: primary-upstream pins the intron start, primary-downstream
    // pins its end. Either way the pinned end fixes the read position, and intronAt reproduces the hinted intron.
    private int junctionPositionFromMate(
            final Candidate candidate, final Side up, final Side down,
            final int upMatchedRead, final boolean primaryIsUpstream)
    {
        for(ChrBaseRegion hint : candidate.mateHintIntrons())
        {
            if(!hint.Chromosome.equals(candidate.chromosome()))
                continue;

            int upLoss = primaryIsUpstream
                    ? up.RefEnd - hint.start() + 1
                    : upMatchedRead - (down.LeadingS + (hint.end() - (down.Start - 1)));

            int readPosition = upMatchedRead - upLoss;
            if(readPosition < down.LeadingS || readPosition > upMatchedRead)
                continue;
            if(up.TrailingM < upLoss || down.LeadingM < readPosition - down.LeadingS)
                continue;

            return readPosition;
        }

        return NO_JUNCTION_POSITION;
    }

    private static final class Side
    {
        int Start;
        List<CigarElement> Cigar;
        int LeadingS;
        int TrailingS;
        int LeadingM;
        int TrailingM;
        int RefEnd;

        private Side(
                final int start, final List<CigarElement> cigar,
                final int leadingS, final int trailingS,
                final int leadingM, final int trailingM, final int refEnd)
        {
            Start = start;
            Cigar = cigar;
            LeadingS = leadingS;
            TrailingS = trailingS;
            LeadingM = leadingM;
            TrailingM = trailingM;
            RefEnd = refEnd;
        }

        static Side of(final int start, final List<CigarElement> cigar)
        {
            return new Side(
                    start, cigar,
                    CigarUtils.leftSoftClipLength(cigar), CigarUtils.rightSoftClipLength(cigar),
                    MergeCigarOps.matchedRun(cigar, false), MergeCigarOps.matchedRun(cigar, true),
                    start + CigarUtils.cigarAlignedLength(cigar) - 1);
        }
    }

    private static final class MergeOutcome
    {
        RejectReason Reject;
        int MergedStart;
        List<CigarElement> MergedCigar;
        ChrBaseRegion IntroducedIntron;
        Supplementary MergedSupp;
        boolean RightExtend;       // which terminal softclip of the primary this merge resolved

        private MergeOutcome(
                final RejectReason reject,
                final int mergedStart, final List<CigarElement> mergedCigar,
                final ChrBaseRegion intron, final Supplementary supp, final boolean rightExtend)
        {
            Reject = reject;
            MergedStart = mergedStart;
            MergedCigar = mergedCigar;
            IntroducedIntron = intron;
            MergedSupp = supp;
            RightExtend = rightExtend;
        }

        static MergeOutcome reject(final RejectReason reason)
        {
            return new MergeOutcome(reason, -1, null, null, null, false);
        }

        static MergeOutcome success(
                final int start, final List<CigarElement> cigar,
                final ChrBaseRegion intron, final Supplementary supp, final boolean rightExtend)
        {
            return new MergeOutcome(null, start, cigar, intron, supp, rightExtend);
        }

        boolean isSuccess()
        {
            return Reject == null;
        }
    }
}
