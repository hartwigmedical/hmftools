package com.hartwig.hmftools.tars.liftback.supplementary;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.liftback.overhang.BoundaryReclaim;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

// Merges a primary's terminal softclip with an annotated-intron-spanning supplementary into a
// single spliced primary (N op). Chains up to SupplementaryConfig.MaxSuppMerges merges per read.
public class SupplementaryResolver
{
    private final AnnotatedJunctionIndex mAnnotatedIndex;
    private final RefSequenceSource mRefSource;
    private final SupplementaryConfig mConfig;
    private final SupplementaryStatistics mStatistics;

    public SupplementaryResolver(final Set<ChrBaseRegion> annotatedJunctions, final SupplementaryConfig config)
    {
        this(new AnnotatedJunctionIndex(annotatedJunctions != null ? annotatedJunctions : new HashSet<>()),
                null, config);
    }

    public SupplementaryResolver(
            final AnnotatedJunctionIndex annotatedIndex, final RefSequenceSource refSource,
            final SupplementaryConfig config)
    {
        mAnnotatedIndex = annotatedIndex != null ? annotatedIndex : new AnnotatedJunctionIndex(new HashSet<>());
        mRefSource = refSource;
        mConfig = config;
        mStatistics = new SupplementaryStatistics(config.MaxSuppMerges);
    }

    public SupplementaryStatistics statistics()
    {
        return mStatistics;
    }

    public SupplementaryResult resolve(final SupplementaryCandidate candidate)
    {
        if(!mConfig.Enabled)
        {
            return SupplementaryResult.noMerge(SupplementaryRejectReason.NO_MATCHING_SUPP);
        }

        if(candidate.supplementaries().isEmpty())
        {
            return tryRefVerifyOnly(candidate);
        }

        mStatistics.countCandidate();

        List<CigarElement> primaryCigar = CigarUtils.cigarElementsFromStr(candidate.primaryCigar());
        if(CigarUtils.hasHardClip(primaryCigar))
        {
            mStatistics.countReject(SupplementaryRejectReason.COMPLEX_CIGAR_SHAPE);
            return SupplementaryResult.noMerge(SupplementaryRejectReason.COMPLEX_CIGAR_SHAPE);
        }

        int primaryStart = candidate.primaryStart();
        List<SupplementaryRecord> remaining = new ArrayList<>(candidate.supplementaries());
        List<Integer> dropped = new ArrayList<>();
        List<ChrBaseRegion> introns = new ArrayList<>();
        SupplementaryRejectReason lastReject = null;
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
                    return SupplementaryResult.noMerge(SupplementaryRejectReason.NO_TERMINAL_SOFTCLIP);
                }
                break;
            }

            MergeOutcome attempt = pickBestSupplementary(
                    candidate, primaryStart, primaryCigar, remaining);

            if(!attempt.isSuccess())
            {
                if(chainDepth == 0)
                {
                    return SupplementaryResult.noMerge(attempt.Reject);
                }
                lastReject = attempt.Reject;
                break;
            }

            primaryStart = attempt.MergedStart;
            primaryCigar = attempt.MergedCigar;
            dropped.add(attempt.MergedSupp.index());
            introns.add(attempt.IntroducedIntron);
            remaining.remove(attempt.MergedSupp);
            ++chainDepth;
        }

        if(chainDepth == 0)
        {
            // Belt and braces: initial-iteration failures should have returned above.
            return SupplementaryResult.noMerge(lastReject != null ? lastReject : SupplementaryRejectReason.NO_MATCHING_SUPP);
        }

        mStatistics.countMergedChain(chainDepth);
        return new SupplementaryResult(
                true, CigarUtils.cigarElementsToStr(primaryCigar), primaryStart,
                dropped, introns, chainDepth, null);
    }

    // No supplementary: look up annotated junctions abutting the softclip boundary and accept if
    // softclipped bases match the candidate exon reference within tolerance.
    private SupplementaryResult tryRefVerifyOnly(final SupplementaryCandidate candidate)
    {
        if(mRefSource == null || candidate.readBases() == null)
        {
            return SupplementaryResult.noMerge(SupplementaryRejectReason.NO_MATCHING_SUPP);
        }

        mStatistics.countCandidate();

        List<CigarElement> primaryCigar = CigarUtils.cigarElementsFromStr(candidate.primaryCigar());
        if(CigarUtils.hasHardClip(primaryCigar))
        {
            mStatistics.countReject(SupplementaryRejectReason.COMPLEX_CIGAR_SHAPE);
            return SupplementaryResult.noMerge(SupplementaryRejectReason.COMPLEX_CIGAR_SHAPE);
        }
        if(primaryCigar.isEmpty())
        {
            return SupplementaryResult.noMerge(SupplementaryRejectReason.NO_TERMINAL_SOFTCLIP);
        }

        boolean trailingS = primaryCigar.get(primaryCigar.size() - 1).getOperator() == CigarOperator.S;
        boolean leadingS = primaryCigar.get(0).getOperator() == CigarOperator.S;
        if(!trailingS && !leadingS)
        {
            return SupplementaryResult.noMerge(SupplementaryRejectReason.NO_TERMINAL_SOFTCLIP);
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
        SupplementaryResult lastFailure = SupplementaryResult.noMerge(SupplementaryRejectReason.REF_VERIFY_MISMATCH_TOO_HIGH);
        for(boolean rightExtend : sides)
        {
            SupplementaryResult sideResult = attemptRefVerifySide(candidate, primaryCigar, rightExtend);
            if(sideResult.merged())
            {
                return sideResult;
            }
            if(sideResult.rejectReason() != SupplementaryRejectReason.REF_VERIFY_NO_CANDIDATE_EXON)
            {
                anyCandidate = true;
                lastFailure = sideResult;
            }
        }

        if(!anyCandidate)
        {
            mStatistics.countReject(SupplementaryRejectReason.REF_VERIFY_NO_CANDIDATE_EXON);
            return SupplementaryResult.noMerge(SupplementaryRejectReason.REF_VERIFY_NO_CANDIDATE_EXON);
        }
        return lastFailure;
    }

    // Ref-verify one terminal softclip against annotated donors/acceptors with boundary snap.
    // Does NOT count the no-candidate reject; tryRefVerifyOnly aggregates across both ends.
    private SupplementaryResult attemptRefVerifySide(
            final SupplementaryCandidate candidate, final List<CigarElement> primaryCigar, final boolean rightExtend)
    {
        if(opAdjacentToSoftClip(primaryCigar, !rightExtend))
        {
            return SupplementaryResult.noMerge(SupplementaryRejectReason.COMPLEX_CIGAR_SHAPE);
        }

        int softclipLen = rightExtend
                ? CigarUtils.rightSoftClipLength(primaryCigar)
                : CigarUtils.leftSoftClipLength(primaryCigar);
        int primaryAnchor = rightExtend
                ? trailingMatchedRun(primaryCigar)
                : leadingMatchedRun(primaryCigar);

        // BWA often over-extends a few bases past the true exon boundary. Snap back up to
        // MaxBoundaryShift (smallest shift first), rolling over-extension into the softclip.
        int primaryRefEnd = candidate.primaryStart() + CigarUtils.cigarAlignedLength(primaryCigar) - 1;
        boolean anyCandidate = false;
        SupplementaryResult lastFailure = SupplementaryResult.noMerge(SupplementaryRejectReason.REF_VERIFY_MISMATCH_TOO_HIGH);

        for(int shift = 0; shift <= mConfig.MaxBoundaryShift; ++shift)
        {
            if(primaryAnchor - shift < 1)
                break;
            List<CigarElement> shiftedCigar = shift == 0
                    ? primaryCigar
                    : shiftBoundaryIntoSoftclip(primaryCigar, shift, rightExtend);
            if(shiftedCigar == null)
                break;

            int boundary = rightExtend ? (primaryRefEnd + 1 - shift) : (candidate.primaryStart() - 1 + shift);
            List<ChrBaseRegion> candidates = rightExtend
                    ? mAnnotatedIndex.introByStart(candidate.chromosome(), boundary)
                    : mAnnotatedIndex.introByEnd(candidate.chromosome(), boundary);
            if(candidates.isEmpty())
                continue;

            anyCandidate = true;
            SupplementaryResult result = verifyAgainstCandidates(
                    candidate, shiftedCigar, candidates, softclipLen + shift, rightExtend);
            if(result.merged())
            {
                return result;
            }
            lastFailure = result;
        }

        if(!anyCandidate)
        {
            return SupplementaryResult.noMerge(SupplementaryRejectReason.REF_VERIFY_NO_CANDIDATE_EXON);
        }
        return lastFailure;
    }

    // Moves 'shift' bases from the exon-edge M into the adjacent softclip to snap back the boundary.
    // Returns null when the edge op isn't M or trimming would leave nothing matched.
    private static List<CigarElement> shiftBoundaryIntoSoftclip(
            final List<CigarElement> cigar, final int shift, final boolean rightExtend)
    {
        List<CigarElement> shifted = new ArrayList<>(cigar);
        int last = shifted.size() - 1;
        int softIdx = rightExtend ? last : 0;
        int matchIdx = rightExtend ? last - 1 : 1;
        if(matchIdx < 0 || matchIdx >= shifted.size())
        {
            return null;
        }

        CigarElement softElement = shifted.get(softIdx);
        CigarElement matchElement = shifted.get(matchIdx);
        if(matchElement.getOperator() != CigarOperator.M && matchElement.getOperator() != CigarOperator.EQ)
        {
            return null;
        }
        if(matchElement.getLength() - shift < 1)
        {
            return null;
        }

        shifted.set(matchIdx, new CigarElement(matchElement.getLength() - shift, matchElement.getOperator()));
        shifted.set(softIdx, new CigarElement(softElement.getLength() + shift, CigarOperator.S));
        return shifted;
    }

    private SupplementaryResult verifyAgainstCandidates(
            final SupplementaryCandidate candidate, final List<CigarElement> primaryCigar,
            final List<ChrBaseRegion> candidates, final int softclipLen, final boolean rightExtend)
    {
        if(candidates.isEmpty())
        {
            mStatistics.countReject(SupplementaryRejectReason.REF_VERIFY_NO_CANDIDATE_EXON);
            return SupplementaryResult.noMerge(SupplementaryRejectReason.REF_VERIFY_NO_CANDIDATE_EXON);
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

            byte[] refBases = mRefSource.getBases(candidate.chromosome(), refStart, refEnd);
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
            mStatistics.countReject(SupplementaryRejectReason.REF_VERIFY_MISMATCH_TOO_HIGH);
            return SupplementaryResult.noMerge(SupplementaryRejectReason.REF_VERIFY_MISMATCH_TOO_HIGH);
        }
        if(ambiguous)
        {
            mStatistics.countReject(SupplementaryRejectReason.REF_VERIFY_AMBIGUOUS);
            return SupplementaryResult.noMerge(SupplementaryRejectReason.REF_VERIFY_AMBIGUOUS);
        }

        return buildRefVerifyMerge(candidate, primaryCigar, chosen, chosenRun, softclipLen, rightExtend);
    }

    private SupplementaryResult buildRefVerifyMerge(
            final SupplementaryCandidate candidate, final List<CigarElement> primaryCigar,
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
        mStatistics.countMergedChain(1);
        return new SupplementaryResult(
                true, CigarUtils.cigarElementsToStr(merged), mergedStart,
                Collections.emptyList(), Collections.singletonList(chosen),
                1, null);
    }

    // Longest run matching ref from the junction-proximal end, scored with the shared bwa-mem model
    // (BoundaryReclaim: match +1, mismatch -4) so ref-verify uses the same matching as collapse/tail-extend.
    // Trailing clips walk forward from index 0 (junction-proximal); leading clips reverse both arrays so the
    // walk runs from the high-index proximal end outward.
    private static int proximalScoringRun(final byte[] softclipBases, final byte[] refBases, final boolean rightExtend)
    {
        if(rightExtend)
        {
            return BoundaryReclaim.maxScoringPrefix(softclipBases, refBases);
        }
        return BoundaryReclaim.maxScoringPrefix(BoundaryReclaim.reversed(softclipBases), BoundaryReclaim.reversed(refBases));
    }

    private static int mismatchesInRun(
            final byte[] softclipBases, final byte[] refBases, final int len, final int run, final boolean proximalAtEnd)
    {
        int mismatches = 0;
        for(int k = 1; k <= run; ++k)
        {
            int index = proximalAtEnd ? (len - k) : (k - 1);
            if(!basesEqualIgnoreCase(softclipBases[index], refBases[index]))
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
        if(mRefSource == null)
        {
            return Tier.NONE;
        }
        byte[] donor = mRefSource.getBases(
                candidateIntron.Chromosome, candidateIntron.start(), candidateIntron.start() + 1);
        byte[] acceptor = mRefSource.getBases(
                candidateIntron.Chromosome, candidateIntron.end() - 1, candidateIntron.end());
        return SpliceMotif.classify(donor, acceptor);
    }

    private static boolean basesEqualIgnoreCase(final byte a, final byte b)
    {
        if(a == b)
        {
            return true;
        }
        return (a & ~0x20) == (b & ~0x20); // mask ASCII case bit; valid for letters, cheap
    }

    private MergeOutcome pickBestSupplementary(
            final SupplementaryCandidate candidate, final int primaryStart,
            final List<CigarElement> primaryCigar, final List<SupplementaryRecord> supps)
    {
        MergeOutcome best = null;
        SupplementaryRejectReason lastReject = null;
        // Ambiguous only when one terminal side draws more than one supp.
        int rightReach = 0;
        int leftReach = 0;

        for(SupplementaryRecord supp : supps)
        {
            MergeOutcome attempt = tryMerge(candidate, primaryStart, primaryCigar, supp);
            if(!attempt.isSuccess())
            {
                mStatistics.countReject(attempt.Reject);
                lastReject = attempt.Reject;
                continue;
            }

            if(attempt.RightExtend)
            {
                ++rightReach;
            }
            else
            {
                ++leftReach;
            }

            if(best == null)
            {
                best = attempt;
            }
            else if(isBetterCandidate(attempt, best))
            {
                best = attempt;
            }
        }

        if(best == null)
        {
            return MergeOutcome.reject(lastReject != null ? lastReject : SupplementaryRejectReason.NO_MATCHING_SUPP);
        }

        if(rightReach > 1 || leftReach > 1)
        {
            mStatistics.countReject(SupplementaryRejectReason.MULTIPLE_SUPPS_IN_REACH);
            return MergeOutcome.reject(SupplementaryRejectReason.MULTIPLE_SUPPS_IN_REACH);
        }

        return best;
    }

    // Higher MAPQ wins, then smaller intron. >1 within reach is rejected by the caller.
    private static boolean isBetterCandidate(final MergeOutcome a, final MergeOutcome b)
    {
        if(a.MergedSupp.mapq() != b.MergedSupp.mapq())
        {
            return a.MergedSupp.mapq() > b.MergedSupp.mapq();
        }
        return intronLength(a.IntroducedIntron) < intronLength(b.IntroducedIntron);
    }

    private static int intronLength(final ChrBaseRegion intron)
    {
        return intron.end() - intron.start() + 1;
    }

    private MergeOutcome tryMerge(
            final SupplementaryCandidate candidate, final int primaryStart,
            final List<CigarElement> primaryCigar, final SupplementaryRecord supp)
    {
        if(!candidate.chromosome().equals(supp.chromosome()))
        {
            return MergeOutcome.reject(SupplementaryRejectReason.DIFFERENT_CHROMOSOME);
        }

        if(candidate.forwardStrand() != supp.forwardStrand())
        {
            return MergeOutcome.reject(SupplementaryRejectReason.OPPOSITE_STRAND);
        }

        List<CigarElement> suppCigar = CigarUtils.cigarElementsFromStr(supp.cigar());
        if(CigarUtils.hasHardClip(suppCigar))
        {
            return MergeOutcome.reject(SupplementaryRejectReason.COMPLEX_CIGAR_SHAPE);
        }

        // ContigTranslator can expand a cross-exon M into M-N-M. If the post-N M coincidentally
        // overlaps the primary's span, clamp the supp to its primary-distal anchor before merge logic.
        int primaryRefEnd = primaryStart + CigarUtils.cigarAlignedLength(primaryCigar) - 1;
        int suppStart = supp.start();
        ClampedSupp clamped = clampSuppToPrimaryBoundary(suppCigar, suppStart, primaryStart, primaryRefEnd);
        if(clamped != null)
        {
            suppCigar = clamped.Cigar;
            suppStart = clamped.Start;
            mStatistics.countSuppClamp();
        }

        int primaryLeadS = CigarUtils.leftSoftClipLength(primaryCigar);
        int primaryTrailS = CigarUtils.rightSoftClipLength(primaryCigar);
        int suppLeadS = CigarUtils.leftSoftClipLength(suppCigar);
        int suppTrailS = CigarUtils.rightSoftClipLength(suppCigar);

        boolean rightExtend = primaryTrailS > 0 && suppLeadS > 0;
        boolean leftExtend = primaryLeadS > 0 && suppTrailS > 0;

        if(!rightExtend && !leftExtend)
        {
            return MergeOutcome.reject(SupplementaryRejectReason.NO_MATCHING_SUPP);
        }

        if(rightExtend && leftExtend)
        {
            int suppRefEnd = suppStart + CigarUtils.cigarAlignedLength(suppCigar) - 1;
            if(suppStart > primaryRefEnd && suppRefEnd >= primaryStart)
            {
                leftExtend = false;
            }
            else if(suppRefEnd < primaryStart && suppStart <= primaryRefEnd)
            {
                rightExtend = false;
            }
            else
            {
                return MergeOutcome.reject(SupplementaryRejectReason.COMPLEX_CIGAR_SHAPE);
            }
        }

        Side primarySide = Side.of(primaryStart, primaryCigar);
        Side suppSide = Side.of(suppStart, suppCigar);
        boolean primaryIsUpstream = rightExtend;
        Side up = primaryIsUpstream ? primarySide : suppSide;
        Side down = primaryIsUpstream ? suppSide : primarySide;

        return mergeJunction(candidate, up, down, primaryIsUpstream, supp);
    }

    // Direction-agnostic merge: validates anchor pair, scores snap points by junction tier, falls
    // back to mate-hint then trust-primary. supp is threaded as result payload only.
    private MergeOutcome mergeJunction(
            final SupplementaryCandidate candidate, final Side up, final Side down,
            final boolean primaryIsUpstream, final SupplementaryRecord supp)
    {
        if(CigarUtils.cigarBaseLength(up.Cigar) != candidate.readLength()
                || CigarUtils.cigarBaseLength(down.Cigar) != candidate.readLength())
            return MergeOutcome.reject(SupplementaryRejectReason.COMPLEX_CIGAR_SHAPE);

        if(opAdjacentToSoftClip(up.Cigar, false) || opAdjacentToSoftClip(down.Cigar, true))
        {
            return MergeOutcome.reject(SupplementaryRejectReason.COMPLEX_CIGAR_SHAPE);
        }

        int upMatchedRead = candidate.readLength() - up.TrailingS;
        int overlap = upMatchedRead - down.LeadingS;

        if(overlap < 0)
        {
            return MergeOutcome.reject(SupplementaryRejectReason.READ_COVERAGE_GAP);
        }
        if(overlap > mConfig.SoftclipTolerance)
        {
            return MergeOutcome.reject(SupplementaryRejectReason.READ_COVERAGE_OVERLAP);
        }
        if(down.Start <= up.RefEnd)
        {
            return MergeOutcome.reject(SupplementaryRejectReason.READ_COVERAGE_OVERLAP);
        }

        // Intron length is invariant under snap point L - check once up front.
        int intronLength = (down.Start - 1 - up.RefEnd) + overlap;
        if(intronLength < mConfig.MinIntronLength)
        {
            return MergeOutcome.reject(SupplementaryRejectReason.INTRON_TOO_SHORT);
        }
        if(intronLength > mConfig.MaxIntronLength)
        {
            return MergeOutcome.reject(SupplementaryRejectReason.INTRON_TOO_LONG);
        }

        // Snap priority: an annotated boundary from the sidecar, else a canonical (then semi-canonical) splice
        // motif; ties within the chosen tier are broken pseudo-randomly but deterministically (seeded by the read),
        // so an ambiguous multi-option junction is distributed yet reproducible. The mate's resolved intron is a
        // strong known-position hint, tried next. With no motif and no hint, the junction is placed at the midpoint
        // of the ambiguous range (rounded down); if even that has no valid anchor split, the merge is not applied and
        // the read keeps whatever placement was already chosen.
        SnapPick pick = scanSnapPoints(candidate, up, down, upMatchedRead, primaryIsUpstream);
        int chosenL = pick.L;
        ChrBaseRegion chosenIntron = pick.Intron;

        if(chosenL == -1)
        {
            SnapPick hinted = scanMateHint(candidate, up, down, upMatchedRead, primaryIsUpstream);
            chosenL = hinted.L;
            chosenIntron = hinted.Intron;
        }
        if(chosenL == -1)
        {
            if(mConfig.AnnotatedOnly)
            {
                return MergeOutcome.reject(SupplementaryRejectReason.NOVEL_JUNCTION);
            }
            // midpoint of the ambiguous read range [down.LeadingS, upMatchedRead], rounded down.
            chosenL = (upMatchedRead + down.LeadingS) / 2;
            int upLossFallback = upMatchedRead - chosenL;
            int downLossFallback = chosenL - down.LeadingS;
            if(up.TrailingM < upLossFallback || down.LeadingM < downLossFallback)
                return MergeOutcome.reject(SupplementaryRejectReason.SHORT_ANCHOR);
            chosenIntron = new ChrBaseRegion(
                    candidate.chromosome(), up.RefEnd - upLossFallback + 1, down.Start + downLossFallback - 1);
        }

        int upLoss = upMatchedRead - chosenL;
        int downLoss = chosenL - down.LeadingS;
        List<CigarElement> merged = buildMergedCigar(up.Cigar, down.Cigar, upLoss, downLoss, intronLength);
        return MergeOutcome.success(up.Start, merged, chosenIntron, supp, primaryIsUpstream);
    }

    // Scores every valid snap point by junction tier and returns one at the highest tier present. When several
    // snap points tie at that tier, the choice is pseudo-random but deterministic (seeded by the read), so an
    // ambiguous junction with multiple equally-good positions is distributed across them yet reproducible run to
    // run. Returns a miss when no snap point reaches a splice motif or annotated boundary.
    private SnapPick scanSnapPoints(
            final SupplementaryCandidate candidate, final Side up, final Side down,
            final int upMatchedRead, final boolean primaryIsUpstream)
    {
        Tier bestTier = Tier.NONE;
        List<Integer> bestLs = new ArrayList<>();
        List<ChrBaseRegion> bestIntrons = new ArrayList<>();

        int startL = primaryIsUpstream ? upMatchedRead : down.LeadingS;
        int endL = primaryIsUpstream ? down.LeadingS : upMatchedRead;
        int step = primaryIsUpstream ? -1 : 1;
        for(int L = startL; primaryIsUpstream ? L >= endL : L <= endL; L += step)
        {
            int upLoss = upMatchedRead - L;
            int downLoss = L - down.LeadingS;
            if(up.TrailingM < upLoss || down.LeadingM < downLoss)
                continue;

            ChrBaseRegion candidateIntron = new ChrBaseRegion(
                    candidate.chromosome(), up.RefEnd - upLoss + 1, down.Start + downLoss - 1);
            Tier tier = classifyJunctionTier(candidateIntron);
            if(tier == Tier.NONE)
                continue;

            int cmp = tier.compareTo(bestTier);
            if(cmp > 0)
            {
                bestTier = tier;
                bestLs.clear();
                bestIntrons.clear();
            }
            if(cmp >= 0)
            {
                bestLs.add(L);
                bestIntrons.add(candidateIntron);
            }
        }

        if(bestLs.isEmpty())
        {
            return SnapPick.miss();
        }
        int index = bestLs.size() == 1 ? 0 : Math.floorMod(snapSeed(candidate, up), bestLs.size());
        return new SnapPick(bestLs.get(index), bestIntrons.get(index));
    }

    // Deterministic per-read seed so an ambiguous multi-option snap is reproducible run to run: hashes the read
    // bases (when available) with the upstream boundary so two reads over the same junction can pick differently.
    private static int snapSeed(final SupplementaryCandidate candidate, final Side up)
    {
        int base = candidate.readBases() != null ? Arrays.hashCode(candidate.readBases()) : 0;
        return 31 * base + up.RefEnd;
    }

    // Use the mate's previously-resolved intron as a hint. Primary-upstream pins intron start;
    // primary-downstream pins intron end.
    private SnapPick scanMateHint(
            final SupplementaryCandidate candidate, final Side up, final Side down,
            final int upMatchedRead, final boolean primaryIsUpstream)
    {
        if(candidate.mateHintIntrons().isEmpty())
        {
            return SnapPick.miss();
        }

        for(ChrBaseRegion hint : candidate.mateHintIntrons())
        {
            if(!hint.Chromosome.equals(candidate.chromosome()))
                continue;
            int upLoss;
            if(primaryIsUpstream)
            {
                upLoss = up.RefEnd - hint.start() + 1;
            }
            else
            {
                upLoss = upMatchedRead - (down.LeadingS + (hint.end() - (down.Start - 1)));
            }
            int L = upMatchedRead - upLoss;
            if(L < down.LeadingS || L > upMatchedRead)
                continue;
            int downLoss = L - down.LeadingS;
            if(up.TrailingM < upLoss || down.LeadingM < downLoss)
                continue;

            int hintedIntronStart = primaryIsUpstream ? hint.start() : (up.RefEnd - upLoss + 1);
            int hintedIntronEnd = primaryIsUpstream ? (down.Start + downLoss - 1) : hint.end();
            return new SnapPick(L, new ChrBaseRegion(candidate.chromosome(), hintedIntronStart, hintedIntronEnd));
        }
        return SnapPick.miss();
    }

    private static List<CigarElement> buildMergedCigar(
            final List<CigarElement> upCigar, final List<CigarElement> downCigar,
            final int upLoss, final int downLoss, final int intronLength)
    {
        List<CigarElement> merged = new ArrayList<>(upCigar.size() + downCigar.size());
        for(int i = 0; i < upCigar.size() - 1; ++i) // upstream ops, excluding trailing S
        {
            if(i == upCigar.size() - 2 && upLoss > 0)
            {
                merged.add(new CigarElement(upCigar.get(i).getLength() - upLoss, upCigar.get(i).getOperator()));
            }
            else
            {
                merged.add(upCigar.get(i));
            }
        }
        merged.add(new CigarElement(intronLength, CigarOperator.N));
        for(int i = 1; i < downCigar.size(); ++i) // downstream ops, excluding leading S
        {
            if(i == 1 && downLoss > 0)
            {
                merged.add(new CigarElement(downCigar.get(i).getLength() - downLoss, downCigar.get(i).getOperator()));
            }
            else
            {
                merged.add(downCigar.get(i));
            }
        }
        return merged;
    }

    // Length of the M/=/X run adjacent to the trailing softclip (e.g. "5M 50S" -> 5).
    private static int trailingMatchedRun(final List<CigarElement> elements)
    {
        for(int i = elements.size() - 1; i >= 0; --i)
        {
            CigarOperator op = elements.get(i).getOperator();
            if(op == CigarOperator.S)
                continue;
            if(op == CigarOperator.M || op == CigarOperator.EQ || op == CigarOperator.X)
            {
                return elements.get(i).getLength();
            }
            return 0;
        }
        return 0;
    }

    private static int leadingMatchedRun(final List<CigarElement> elements)
    {
        for(int i = 0; i < elements.size(); ++i)
        {
            CigarOperator op = elements.get(i).getOperator();
            if(op == CigarOperator.S)
                continue;
            if(op == CigarOperator.M || op == CigarOperator.EQ || op == CigarOperator.X)
            {
                return elements.get(i).getLength();
            }
            return 0;
        }
        return 0;
    }

    private static boolean opAdjacentToSoftClip(final List<CigarElement> elements, final boolean leadingSide)
    {
        if(elements.size() < 2)
        {
            return false;
        }
        if(leadingSide)
        {
            if(elements.get(0).getOperator() != CigarOperator.S)
            {
                return false;
            }
            CigarOperator next = elements.get(1).getOperator();
            return next == CigarOperator.I || next == CigarOperator.D;
        }
        else
        {
            int last = elements.size() - 1;
            if(elements.get(last).getOperator() != CigarOperator.S)
            {
                return false;
            }
            CigarOperator prev = elements.get(last - 1).getOperator();
            return prev == CigarOperator.I || prev == CigarOperator.D;
        }
    }

    // Replaces the supp's post-N tail (or pre-N head) with a softclip when the M-after-N overlaps
    // the primary's span (ContigTranslator M-N-M artifact double-counting read bases). Returns null
    // when no internal N or no overlap.
    private static ClampedSupp clampSuppToPrimaryBoundary(
            final List<CigarElement> suppCigar, final int suppStart,
            final int primaryStart, final int primaryRefEnd)
    {
        boolean hasInternalN = false;
        for(CigarElement e : suppCigar)
        {
            if(e.getOperator() == CigarOperator.N)
            {
                hasInternalN = true;
                break;
            }
        }
        if(!hasInternalN)
        {
            return null;
        }

        int suppRefEnd = suppStart + CigarUtils.cigarAlignedLength(suppCigar) - 1;
        if(suppStart < primaryStart)
        {
            return clampUpstreamSupp(suppCigar, suppStart, primaryStart);
        }
        if(suppRefEnd > primaryRefEnd)
        {
            return clampDownstreamSupp(suppCigar, suppStart, primaryRefEnd);
        }
        return null;
    }

    private static ClampedSupp clampUpstreamSupp(
            final List<CigarElement> suppCigar, final int suppStart, final int primaryStart)
    {
        int refCursor = suppStart;
        int readCursor = 0;
        int firstNIdx = -1;
        int readAtFirstN = 0;

        for(int i = 0; i < suppCigar.size(); i++)
        {
            CigarElement e = suppCigar.get(i);
            CigarOperator op = e.getOperator();
            int len = e.getLength();

            if(op == CigarOperator.N && firstNIdx == -1)
            {
                firstNIdx = i;
                readAtFirstN = readCursor;
            }

            if(firstNIdx != -1 && refCursor >= primaryStart
                    && (op == CigarOperator.M || op == CigarOperator.EQ || op == CigarOperator.X))
            {
                return cutTailAt(suppCigar, suppStart, firstNIdx, readAtFirstN);
            }

            if(op.consumesReferenceBases())
            {
                refCursor += len;
            }
            if(op.consumesReadBases())
            {
                readCursor += len;
            }
        }
        return null;
    }

    private static ClampedSupp clampDownstreamSupp(
            final List<CigarElement> suppCigar, final int suppStart, final int primaryRefEnd)
    {
        int refCursor = suppStart;
        int readCursor = 0;
        int lastBadNIdx = -1;
        int readAfterLastBadN = 0;
        int refAfterLastBadN = -1;
        boolean prevMInsidePrimary = false;

        for(int i = 0; i < suppCigar.size(); i++)
        {
            CigarElement e = suppCigar.get(i);
            CigarOperator op = e.getOperator();
            int len = e.getLength();

            if(op == CigarOperator.M || op == CigarOperator.EQ || op == CigarOperator.X)
            {
                int mRefEnd = refCursor + len - 1;
                prevMInsidePrimary = mRefEnd <= primaryRefEnd;
            }

            if(op == CigarOperator.N && prevMInsidePrimary)
            {
                lastBadNIdx = i;
                refAfterLastBadN = refCursor + len; // just after this N op
                readAfterLastBadN = readCursor;
            }

            if(op.consumesReferenceBases())
            {
                refCursor += len;
            }
            if(op.consumesReadBases())
            {
                readCursor += len;
            }
        }

        if(lastBadNIdx == -1 || refAfterLastBadN < 0)
        {
            return null;
        }
        return cutHeadAt(suppCigar, refAfterLastBadN, lastBadNIdx, readAfterLastBadN);
    }

    private static ClampedSupp cutTailAt(
            final List<CigarElement> suppCigar, final int suppStart,
            final int cutIdx, final int readBeforeCut)
    {
        List<CigarElement> trimmed = new ArrayList<>(cutIdx + 1);
        for(int i = 0; i < cutIdx; i++)
        {
            trimmed.add(suppCigar.get(i));
        }
        int totalRead = CigarUtils.cigarBaseLength(suppCigar);
        int trailingS = totalRead - readBeforeCut;
        if(trailingS > 0)
        {
            trimmed.add(new CigarElement(trailingS, CigarOperator.S));
        }
        return new ClampedSupp(suppStart, trimmed);
    }

    private static ClampedSupp cutHeadAt(
            final List<CigarElement> suppCigar, final int newStart,
            final int cutIdx, final int readBeforeCut)
    {
        List<CigarElement> trimmed = new ArrayList<>(suppCigar.size() - cutIdx);
        if(readBeforeCut > 0)
        {
            trimmed.add(new CigarElement(readBeforeCut, CigarOperator.S));
        }
        for(int i = cutIdx + 1; i < suppCigar.size(); i++)
        {
            trimmed.add(suppCigar.get(i));
        }
        return new ClampedSupp(newStart, trimmed);
    }

    private static final class ClampedSupp
    {
        int Start;
        List<CigarElement> Cigar;

        ClampedSupp(final int start, final List<CigarElement> cigar)
        {
            Start = start;
            Cigar = cigar;
        }
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
            return new Side(start, cigar,
                    CigarUtils.leftSoftClipLength(cigar), CigarUtils.rightSoftClipLength(cigar),
                    leadingMatchedRun(cigar), trailingMatchedRun(cigar),
                    start + CigarUtils.cigarAlignedLength(cigar) - 1);
        }
    }

    private static final class SnapPick
    {
        int L;
        ChrBaseRegion Intron;

        SnapPick(final int snapPoint, final ChrBaseRegion intron)
        {
            L = snapPoint;
            Intron = intron;
        }

        static SnapPick miss()
        {
            return new SnapPick(-1, null);
        }
    }

    private static final class MergeOutcome
    {
        SupplementaryRejectReason Reject;
        int MergedStart;
        List<CigarElement> MergedCigar;
        ChrBaseRegion IntroducedIntron;
        SupplementaryRecord MergedSupp;
        boolean RightExtend;       // which terminal softclip of the primary this merge resolved

        private MergeOutcome(
                final SupplementaryRejectReason reject,
                final int mergedStart, final List<CigarElement> mergedCigar,
                final ChrBaseRegion intron, final SupplementaryRecord supp, final boolean rightExtend)
        {
            Reject = reject;
            MergedStart = mergedStart;
            MergedCigar = mergedCigar;
            IntroducedIntron = intron;
            MergedSupp = supp;
            RightExtend = rightExtend;
        }

        static MergeOutcome reject(final SupplementaryRejectReason reason)
        {
            return new MergeOutcome(reason, -1, null, null, null, false);
        }

        static MergeOutcome success(
                final int start, final List<CigarElement> cigar,
                final ChrBaseRegion intron, final SupplementaryRecord supp, final boolean rightExtend)
        {
            return new MergeOutcome(null, start, cigar, intron, supp, rightExtend);
        }

        boolean isSuccess()
        {
            return Reject == null;
        }
    }
}
