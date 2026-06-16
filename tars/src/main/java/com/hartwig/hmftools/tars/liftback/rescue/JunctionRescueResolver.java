package com.hartwig.hmftools.tars.liftback.rescue;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.common.SpliceCommon;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

// Merges a primary's terminal softclip with an annotated-intron-spanning supplementary into a
// single spliced primary (N op). Chains up to RescueConfig.MaxChainDepth merges per read.
public class JunctionRescueResolver
{
    private final AnnotatedJunctionIndex mAnnotatedIndex;
    private final RefSequenceSource mRefSource;
    private final RescueConfig mConfig;
    private final RescueStatistics mStatistics;

    public JunctionRescueResolver(final Set<ChrBaseRegion> annotatedJunctions, final RescueConfig config)
    {
        this(new AnnotatedJunctionIndex(annotatedJunctions != null ? annotatedJunctions : new HashSet<>()),
                null, config);
    }

    public JunctionRescueResolver(
            final AnnotatedJunctionIndex annotatedIndex, final RefSequenceSource refSource,
            final RescueConfig config)
    {
        mAnnotatedIndex = annotatedIndex != null ? annotatedIndex : new AnnotatedJunctionIndex(new HashSet<>());
        mRefSource = refSource;
        mConfig = config;
        mStatistics = new RescueStatistics(config.MaxChainDepth);
    }

    public RescueStatistics statistics()
    {
        return mStatistics;
    }

    public RescueResult resolve(final RescueCandidate candidate)
    {
        if(!mConfig.Enabled)
            return RescueResult.noMerge(RescueRejectReason.NO_MATCHING_SUPP);

        if(candidate.supplementaries().isEmpty())
            return tryRefVerifyOnly(candidate);

        mStatistics.countCandidate();

        // Low-confidence placement isn't a trustworthy anchor for building a spliced read.
        if(candidate.primaryMapq() < mConfig.MinPrimaryMapq)
        {
            mStatistics.countReject(RescueRejectReason.LOW_PRIMARY_MAPQ);
            return RescueResult.noMerge(RescueRejectReason.LOW_PRIMARY_MAPQ);
        }

        List<CigarElement> primaryCigar = CigarUtils.cigarElementsFromStr(candidate.primaryCigar());
        if(CigarUtils.hasHardClip(primaryCigar))
        {
            mStatistics.countReject(RescueRejectReason.COMPLEX_CIGAR_SHAPE);
            return RescueResult.noMerge(RescueRejectReason.COMPLEX_CIGAR_SHAPE);
        }

        // Fold tx-contig lift artifact micro-junctions into the softclip so the merge anchors on the
        // real matched block. Only affects the working cigar; emitted cigar changes only on accept.
        primaryCigar = foldFabricatedTerminalMicroJunctions(primaryCigar, SpliceCommon.MIN_JUNCTION_ANCHOR);

        int primaryStart = candidate.primaryStart();
        final List<RescueSupplementary> remaining = new ArrayList<>(candidate.supplementaries());
        final List<Integer> dropped = new ArrayList<>();
        final List<ChrBaseRegion> introns = new ArrayList<>();
        RescueRejectReason lastReject = null;
        int chainDepth = 0;

        while(chainDepth < mConfig.MaxChainDepth && !remaining.isEmpty())
        {
            final boolean primaryHasLeadingS = !primaryCigar.isEmpty()
                    && primaryCigar.get(0).getOperator() == CigarOperator.S;
            final boolean primaryHasTrailingS = !primaryCigar.isEmpty()
                    && primaryCigar.get(primaryCigar.size() - 1).getOperator() == CigarOperator.S;
            if(!primaryHasLeadingS && !primaryHasTrailingS)
            {
                if(chainDepth == 0)
                    return RescueResult.noMerge(RescueRejectReason.NO_TERMINAL_SOFTCLIP);
                break;
            }

            final MergeOutcome attempt = pickBestSupplementary(
                    candidate, primaryStart, primaryCigar, remaining);

            if(!attempt.isSuccess())
            {
                if(chainDepth == 0)
                    return RescueResult.noMerge(attempt.Reject);
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
            return RescueResult.noMerge(lastReject != null ? lastReject : RescueRejectReason.NO_MATCHING_SUPP);
        }

        mStatistics.countMergedChain(chainDepth);
        return new RescueResult(
                true, CigarUtils.cigarElementsToStr(primaryCigar), primaryStart,
                dropped, introns, chainDepth, null);
    }

    // No supplementary: look up annotated junctions abutting the softclip boundary and accept if
    // softclipped bases match the candidate exon reference within tolerance.
    private RescueResult tryRefVerifyOnly(final RescueCandidate candidate)
    {
        if(mRefSource == null || candidate.readBases() == null)
            return RescueResult.noMerge(RescueRejectReason.NO_MATCHING_SUPP);

        mStatistics.countCandidate();

        final List<CigarElement> primaryCigar = CigarUtils.cigarElementsFromStr(candidate.primaryCigar());
        if(CigarUtils.hasHardClip(primaryCigar))
        {
            mStatistics.countReject(RescueRejectReason.COMPLEX_CIGAR_SHAPE);
            return RescueResult.noMerge(RescueRejectReason.COMPLEX_CIGAR_SHAPE);
        }
        if(primaryCigar.isEmpty())
            return RescueResult.noMerge(RescueRejectReason.NO_TERMINAL_SOFTCLIP);

        final boolean trailingS = primaryCigar.get(primaryCigar.size() - 1).getOperator() == CigarOperator.S;
        final boolean leadingS = primaryCigar.get(0).getOperator() == CigarOperator.S;
        if(!trailingS && !leadingS)
            return RescueResult.noMerge(RescueRejectReason.NO_TERMINAL_SOFTCLIP);

        // Both-end clips: try longer clip first (more likely the splice tail). First win is taken.
        final boolean[] sides;
        if(trailingS && leadingS)
        {
            final boolean trailingFirst =
                    CigarUtils.rightSoftClipLength(primaryCigar) >= CigarUtils.leftSoftClipLength(primaryCigar);
            sides = new boolean[] { trailingFirst, !trailingFirst };
        }
        else
        {
            sides = new boolean[] { trailingS };
        }

        boolean anyCandidate = false;
        RescueResult lastFailure = RescueResult.noMerge(RescueRejectReason.REF_VERIFY_MISMATCH_TOO_HIGH);
        for(boolean rightExtend : sides)
        {
            final RescueResult sideResult = attemptRefVerifySide(candidate, primaryCigar, rightExtend);
            if(sideResult.merged())
                return sideResult;
            if(sideResult.rejectReason() != RescueRejectReason.REF_VERIFY_NO_CANDIDATE_EXON)
            {
                anyCandidate = true;
                lastFailure = sideResult;
            }
        }

        if(!anyCandidate)
        {
            mStatistics.countReject(RescueRejectReason.REF_VERIFY_NO_CANDIDATE_EXON);
            return RescueResult.noMerge(RescueRejectReason.REF_VERIFY_NO_CANDIDATE_EXON);
        }
        return lastFailure;
    }

    // Ref-verify one terminal softclip against annotated donors/acceptors with boundary snap.
    // Does NOT count the no-candidate reject; tryRefVerifyOnly aggregates across both ends.
    private RescueResult attemptRefVerifySide(
            final RescueCandidate candidate, final List<CigarElement> primaryCigar, final boolean rightExtend)
    {
        if(opAdjacentToSoftClip(primaryCigar, !rightExtend))
            return RescueResult.noMerge(RescueRejectReason.COMPLEX_CIGAR_SHAPE);

        final int softclipLen = rightExtend
                ? CigarUtils.rightSoftClipLength(primaryCigar)
                : CigarUtils.leftSoftClipLength(primaryCigar);
        if(softclipLen < mConfig.MinAnchorOverhang)
            return RescueResult.noMerge(RescueRejectReason.SHORT_ANCHOR);
        final int primaryAnchor = rightExtend
                ? trailingMatchedRun(primaryCigar)
                : leadingMatchedRun(primaryCigar);
        if(primaryAnchor < mConfig.MinAnchorOverhang)
            return RescueResult.noMerge(RescueRejectReason.SHORT_ANCHOR);

        // BWA often over-extends a few bases past the true exon boundary. Snap back up to
        // MaxBoundaryShift (smallest shift first), rolling over-extension into the softclip.
        final int primaryRefEnd = candidate.primaryStart() + CigarUtils.cigarAlignedLength(primaryCigar) - 1;
        boolean anyCandidate = false;
        RescueResult lastFailure = RescueResult.noMerge(RescueRejectReason.REF_VERIFY_MISMATCH_TOO_HIGH);

        for(int shift = 0; shift <= mConfig.MaxBoundaryShift; ++shift)
        {
            if(primaryAnchor - shift < mConfig.MinAnchorOverhang)
                break;
            final List<CigarElement> shiftedCigar = shift == 0
                    ? primaryCigar
                    : shiftBoundaryIntoSoftclip(primaryCigar, shift, rightExtend);
            if(shiftedCigar == null)
                break;

            final int boundary = rightExtend ? (primaryRefEnd + 1 - shift) : (candidate.primaryStart() - 1 + shift);
            final List<ChrBaseRegion> candidates = rightExtend
                    ? mAnnotatedIndex.introByStart(candidate.chromosome(), boundary)
                    : mAnnotatedIndex.introByEnd(candidate.chromosome(), boundary);
            if(candidates.isEmpty())
                continue;

            anyCandidate = true;
            final RescueResult result = verifyAgainstCandidates(
                    candidate, shiftedCigar, candidates, softclipLen + shift, rightExtend);
            if(result.merged())
                return result;
            lastFailure = result;
        }

        if(!anyCandidate)
            return RescueResult.noMerge(RescueRejectReason.REF_VERIFY_NO_CANDIDATE_EXON);
        return lastFailure;
    }

    // Moves 'shift' bases from the exon-edge M into the adjacent softclip to snap back the boundary.
    // Returns null when the edge op isn't M or trimming would leave nothing matched.
    private static List<CigarElement> shiftBoundaryIntoSoftclip(
            final List<CigarElement> cigar, final int shift, final boolean rightExtend)
    {
        final List<CigarElement> shifted = new ArrayList<>(cigar);
        final int last = shifted.size() - 1;
        final int softIdx = rightExtend ? last : 0;
        final int matchIdx = rightExtend ? last - 1 : 1;
        if(matchIdx < 0 || matchIdx >= shifted.size())
            return null;

        final CigarElement softElement = shifted.get(softIdx);
        final CigarElement matchElement = shifted.get(matchIdx);
        if(matchElement.getOperator() != CigarOperator.M && matchElement.getOperator() != CigarOperator.EQ)
            return null;
        if(matchElement.getLength() - shift < 1)
            return null;

        shifted.set(matchIdx, new CigarElement(matchElement.getLength() - shift, matchElement.getOperator()));
        shifted.set(softIdx, new CigarElement(softElement.getLength() + shift, CigarOperator.S));
        return shifted;
    }

    private RescueResult verifyAgainstCandidates(
            final RescueCandidate candidate, final List<CigarElement> primaryCigar,
            final List<ChrBaseRegion> candidates, final int softclipLen, final boolean rightExtend)
    {
        if(candidates.isEmpty())
        {
            mStatistics.countReject(RescueRejectReason.REF_VERIFY_NO_CANDIDATE_EXON);
            return RescueResult.noMerge(RescueRejectReason.REF_VERIFY_NO_CANDIDATE_EXON);
        }

        final byte[] readBases = candidate.readBases();
        final byte[] softclipBases = new byte[softclipLen];
        if(rightExtend)
            System.arraycopy(readBases, readBases.length - softclipLen, softclipBases, 0, softclipLen);
        else
            System.arraycopy(readBases, 0, softclipBases, 0, softclipLen);

        // Match from the junction-proximal end (10% mismatch floor). The outer softclip residual stays
        // soft-clipped since bwa often carries adapter/low-quality bases there (e.g. 19S132M ->
        // 6S15M..N..130M, not 21M..N..130M). Proximal end is the high index for a leading softclip.
        final boolean proximalAtEnd = !rightExtend;
        ChrBaseRegion chosen = null;
        int chosenRun = 0;
        int chosenMismatches = Integer.MAX_VALUE;
        boolean ambiguous = false;

        for(ChrBaseRegion candidateIntron : candidates)
        {
            final int intronLength = candidateIntron.end() - candidateIntron.start() + 1;
            if(intronLength < mConfig.MinIntronLength || intronLength > mConfig.MaxIntronLength)
                continue;

            final int refStart;
            final int refEnd;
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

            final byte[] refBases = mRefSource.getBases(candidate.chromosome(), refStart, refEnd);
            if(refBases == null || refBases.length != softclipLen)
                continue;

            final int run = longestProximalMatchRun(softclipBases, refBases, softclipLen, proximalAtEnd);
            if(run == 0)
                continue;
            final int mismatches = mismatchesInRun(softclipBases, refBases, softclipLen, run, proximalAtEnd);

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
            mStatistics.countReject(RescueRejectReason.REF_VERIFY_MISMATCH_TOO_HIGH);
            return RescueResult.noMerge(RescueRejectReason.REF_VERIFY_MISMATCH_TOO_HIGH);
        }
        if(ambiguous)
        {
            mStatistics.countReject(RescueRejectReason.REF_VERIFY_AMBIGUOUS);
            return RescueResult.noMerge(RescueRejectReason.REF_VERIFY_AMBIGUOUS);
        }
        // Partial match (outer residual stays clipped) requires a longer proximal anchor than a full match.
        final boolean fullMatch = chosenRun == softclipLen;
        if(!fullMatch && chosenRun < mConfig.MinPartialMatchRun)
        {
            mStatistics.countReject(RescueRejectReason.REF_VERIFY_SHORT_PARTIAL_RUN);
            return RescueResult.noMerge(RescueRejectReason.REF_VERIFY_SHORT_PARTIAL_RUN);
        }

        return buildRefVerifyMerge(candidate, primaryCigar, chosen, chosenRun, softclipLen, rightExtend);
    }

    private RescueResult buildRefVerifyMerge(
            final RescueCandidate candidate, final List<CigarElement> primaryCigar,
            final ChrBaseRegion chosen, final int matchedRun, final int softclipLen, final boolean rightExtend)
    {
        final int intronLength = chosen.end() - chosen.start() + 1;
        final int residualSoftclip = softclipLen - matchedRun;
        final List<CigarElement> merged = new ArrayList<>(primaryCigar.size() + 3);
        if(rightExtend)
        {
            for(int i = 0; i < primaryCigar.size() - 1; ++i)
                merged.add(primaryCigar.get(i));
            merged.add(new CigarElement(intronLength, CigarOperator.N));
            merged.add(new CigarElement(matchedRun, CigarOperator.M));
            if(residualSoftclip > 0)
                merged.add(new CigarElement(residualSoftclip, CigarOperator.S));
        }
        else
        {
            if(residualSoftclip > 0)
                merged.add(new CigarElement(residualSoftclip, CigarOperator.S));
            merged.add(new CigarElement(matchedRun, CigarOperator.M));
            merged.add(new CigarElement(intronLength, CigarOperator.N));
            for(int i = 1; i < primaryCigar.size(); ++i)
                merged.add(primaryCigar.get(i));
        }

        final int mergedStart = rightExtend ? candidate.primaryStart() : (chosen.start() - matchedRun);
        mStatistics.countMergedChain(1);
        return new RescueResult(
                true, CigarUtils.cigarElementsToStr(merged), mergedStart,
                Collections.emptyList(), Collections.singletonList(chosen),
                1, null);
    }

    // Longest run matching ref from the junction-proximal end, within a cumulative 10% mismatch floor.
    // proximalAtEnd=true walks inward from the high index; ref[i] pairs with softclip[i].
    private static int longestProximalMatchRun(
            final byte[] softclipBases, final byte[] refBases, final int len, final boolean proximalAtEnd)
    {
        int mismatches = 0;
        int best = 0;
        for(int k = 1; k <= len; ++k)
        {
            final int idx = proximalAtEnd ? (len - k) : (k - 1);
            final boolean match = basesEqualIgnoreCase(softclipBases[idx], refBases[idx]);
            if(!match)
                ++mismatches;
            if(mismatches > k / 10)
                break;
            if(match)
                best = k;
        }
        return best;
    }

    private static int mismatchesInRun(
            final byte[] softclipBases, final byte[] refBases, final int len, final int run, final boolean proximalAtEnd)
    {
        int mismatches = 0;
        for(int k = 1; k <= run; ++k)
        {
            final int idx = proximalAtEnd ? (len - k) : (k - 1);
            if(!basesEqualIgnoreCase(softclipBases[idx], refBases[idx]))
                ++mismatches;
        }
        return mismatches;
    }

    private int classifyJunctionTier(final ChrBaseRegion candidateIntron)
    {
        if(mAnnotatedIndex.contains(candidateIntron))
            return SpliceMotif.TIER_ANNOTATED;
        if(mRefSource == null)
            return SpliceMotif.TIER_NONE;
        final byte[] donor = mRefSource.getBases(
                candidateIntron.Chromosome, candidateIntron.start(), candidateIntron.start() + 1);
        final byte[] acceptor = mRefSource.getBases(
                candidateIntron.Chromosome, candidateIntron.end() - 1, candidateIntron.end());
        return SpliceMotif.classify(donor, acceptor);
    }

    private static boolean basesEqualIgnoreCase(final byte a, final byte b)
    {
        if(a == b)
            return true;
        return (a & ~0x20) == (b & ~0x20); // mask ASCII case bit; valid for letters, cheap
    }

    private MergeOutcome pickBestSupplementary(
            final RescueCandidate candidate, final int primaryStart,
            final List<CigarElement> primaryCigar, final List<RescueSupplementary> supps)
    {
        MergeOutcome best = null;
        RescueRejectReason lastReject = null;
        // Ambiguous only when one terminal side draws more than one supp.
        int rightReach = 0;
        int leftReach = 0;

        for(RescueSupplementary supp : supps)
        {
            final MergeOutcome attempt = tryMerge(candidate, primaryStart, primaryCigar, supp);
            if(!attempt.isSuccess())
            {
                mStatistics.countReject(attempt.Reject);
                lastReject = attempt.Reject;
                continue;
            }

            if(attempt.RightExtend)
                ++rightReach;
            else
                ++leftReach;

            if(best == null)
                best = attempt;
            else if(isBetterCandidate(attempt, best))
                best = attempt;
        }

        if(best == null)
            return MergeOutcome.reject(lastReject != null ? lastReject : RescueRejectReason.NO_MATCHING_SUPP);

        if(rightReach > 1 || leftReach > 1)
        {
            mStatistics.countReject(RescueRejectReason.MULTIPLE_SUPPS_IN_REACH);
            return MergeOutcome.reject(RescueRejectReason.MULTIPLE_SUPPS_IN_REACH);
        }

        return best;
    }

    // Higher MAPQ wins, then smaller intron. >1 within reach is rejected by the caller.
    private static boolean isBetterCandidate(final MergeOutcome a, final MergeOutcome b)
    {
        if(a.MergedSupp.mapq() != b.MergedSupp.mapq())
            return a.MergedSupp.mapq() > b.MergedSupp.mapq();
        return intronLength(a.IntroducedIntron) < intronLength(b.IntroducedIntron);
    }

    private static int intronLength(final ChrBaseRegion intron)
    {
        return intron.end() - intron.start() + 1;
    }

    private MergeOutcome tryMerge(
            final RescueCandidate candidate, final int primaryStart,
            final List<CigarElement> primaryCigar, final RescueSupplementary supp)
    {
        if(!candidate.chromosome().equals(supp.chromosome()))
            return MergeOutcome.reject(RescueRejectReason.DIFFERENT_CHROMOSOME);

        if(candidate.forwardStrand() != supp.forwardStrand())
            return MergeOutcome.reject(RescueRejectReason.OPPOSITE_STRAND);

        List<CigarElement> suppCigar = CigarUtils.cigarElementsFromStr(supp.cigar());
        if(CigarUtils.hasHardClip(suppCigar))
            return MergeOutcome.reject(RescueRejectReason.COMPLEX_CIGAR_SHAPE);

        // ContigTranslator can expand a cross-exon M into M-N-M. If the post-N M coincidentally
        // overlaps the primary's span, clamp the supp to its primary-distal anchor before merge logic.
        final int primaryRefEnd = primaryStart + CigarUtils.cigarAlignedLength(primaryCigar) - 1;
        int suppStart = supp.start();
        final ClampedSupp clamped = clampSuppToPrimaryBoundary(suppCigar, suppStart, primaryStart, primaryRefEnd);
        if(clamped != null)
        {
            suppCigar = clamped.Cigar;
            suppStart = clamped.Start;
            mStatistics.countSuppClamp();
        }

        final int primaryLeadS = CigarUtils.leftSoftClipLength(primaryCigar);
        final int primaryTrailS = CigarUtils.rightSoftClipLength(primaryCigar);
        final int suppLeadS = CigarUtils.leftSoftClipLength(suppCigar);
        final int suppTrailS = CigarUtils.rightSoftClipLength(suppCigar);

        boolean rightExtend = primaryTrailS > 0 && suppLeadS > 0;
        boolean leftExtend = primaryLeadS > 0 && suppTrailS > 0;

        if(!rightExtend && !leftExtend)
            return MergeOutcome.reject(RescueRejectReason.NO_MATCHING_SUPP);

        if(rightExtend && leftExtend)
        {
            final int suppRefEnd = suppStart + CigarUtils.cigarAlignedLength(suppCigar) - 1;
            if(suppStart > primaryRefEnd && suppRefEnd >= primaryStart)
                leftExtend = false;
            else if(suppRefEnd < primaryStart && suppStart <= primaryRefEnd)
                rightExtend = false;
            else
                return MergeOutcome.reject(RescueRejectReason.COMPLEX_CIGAR_SHAPE);
        }

        final Side primarySide = Side.of(primaryStart, primaryCigar);
        final Side suppSide = Side.of(suppStart, suppCigar);
        final boolean primaryIsUpstream = rightExtend;
        final Side up = primaryIsUpstream ? primarySide : suppSide;
        final Side down = primaryIsUpstream ? suppSide : primarySide;

        return mergeJunction(candidate, up, down, primaryIsUpstream, supp);
    }

    // Direction-agnostic merge: validates anchor pair, scores snap points by junction tier, falls
    // back to mate-hint then trust-primary. supp is threaded as result payload only.
    private MergeOutcome mergeJunction(
            final RescueCandidate candidate, final Side up, final Side down,
            final boolean primaryIsUpstream, final RescueSupplementary supp)
    {
        if(CigarUtils.cigarBaseLength(up.Cigar) != candidate.readLength()
                || CigarUtils.cigarBaseLength(down.Cigar) != candidate.readLength())
            return MergeOutcome.reject(RescueRejectReason.COMPLEX_CIGAR_SHAPE);

        if(opAdjacentToSoftClip(up.Cigar, false) || opAdjacentToSoftClip(down.Cigar, true))
            return MergeOutcome.reject(RescueRejectReason.COMPLEX_CIGAR_SHAPE);

        final int upMatchedRead = candidate.readLength() - up.TrailingS;
        final int overlap = upMatchedRead - down.LeadingS;

        if(overlap < 0)
            return MergeOutcome.reject(RescueRejectReason.READ_COVERAGE_GAP);
        if(overlap > mConfig.SoftclipTolerance)
            return MergeOutcome.reject(RescueRejectReason.READ_COVERAGE_OVERLAP);
        if(down.Start <= up.RefEnd)
            return MergeOutcome.reject(RescueRejectReason.READ_COVERAGE_OVERLAP);

        // Intron length is invariant under snap point L - check once up front.
        final int intronLength = (down.Start - 1 - up.RefEnd) + overlap;
        if(intronLength < mConfig.MinIntronLength)
            return MergeOutcome.reject(RescueRejectReason.INTRON_TOO_SHORT);
        if(intronLength > mConfig.MaxIntronLength)
            return MergeOutcome.reject(RescueRejectReason.INTRON_TOO_LONG);

        if(up.TrailingM < mConfig.MinAnchorOverhang || down.LeadingM < mConfig.MinAnchorOverhang)
            return MergeOutcome.reject(RescueRejectReason.SHORT_ANCHOR);

        // Highest-tier snap wins; within a tier, largest min(upAnchor, downAnchor). Trust-primary
        // is visited first so it wins ties (upLoss=0 / downLoss=0 depending on direction).
        final SnapPick pick = scanSnapPoints(candidate, up, down, upMatchedRead, primaryIsUpstream);
        int chosenL = pick.L;
        ChrBaseRegion chosenIntron = pick.Intron;

        if(chosenL == -1)
        {
            final SnapPick hinted = scanMateHint(candidate, up, down, upMatchedRead, primaryIsUpstream);
            chosenL = hinted.L;
            chosenIntron = hinted.Intron;
        }
        if(chosenL == -1)
        {
            if(mConfig.AnnotatedOnly)
                return MergeOutcome.reject(RescueRejectReason.NOVEL_JUNCTION);
            chosenL = primaryIsUpstream ? upMatchedRead : down.LeadingS;
            final int upLossFb = upMatchedRead - chosenL;
            final int downLossFb = chosenL - down.LeadingS;
            if(up.TrailingM - upLossFb < mConfig.MinAnchorOverhang || up.TrailingM < upLossFb
                    || down.LeadingM - downLossFb < mConfig.MinAnchorOverhang || down.LeadingM < downLossFb)
                return MergeOutcome.reject(RescueRejectReason.SHORT_ANCHOR);
            chosenIntron = new ChrBaseRegion(candidate.chromosome(), up.RefEnd + 1, down.Start - 1);
        }

        final int upLoss = upMatchedRead - chosenL;
        final int downLoss = chosenL - down.LeadingS;
        final List<CigarElement> merged = buildMergedCigar(up.Cigar, down.Cigar, upLoss, downLoss, intronLength);
        return MergeOutcome.success(up.Start, merged, chosenIntron, supp, primaryIsUpstream);
    }

    private SnapPick scanSnapPoints(
            final RescueCandidate candidate, final Side up, final Side down,
            final int upMatchedRead, final boolean primaryIsUpstream)
    {
        int chosenL = -1;
        ChrBaseRegion chosenIntron = null;
        int chosenTier = SpliceMotif.TIER_NONE;
        int chosenMinAnchor = -1;

        final int startL = primaryIsUpstream ? upMatchedRead : down.LeadingS;
        final int endL = primaryIsUpstream ? down.LeadingS : upMatchedRead;
        final int step = primaryIsUpstream ? -1 : 1;
        for(int L = startL; primaryIsUpstream ? L >= endL : L <= endL; L += step)
        {
            final int upLoss = upMatchedRead - L;
            final int downLoss = L - down.LeadingS;
            if(up.TrailingM - upLoss < mConfig.MinAnchorOverhang)
                continue;
            if(down.LeadingM - downLoss < mConfig.MinAnchorOverhang)
                continue;
            if(up.TrailingM < upLoss || down.LeadingM < downLoss)
                continue;

            final ChrBaseRegion candidateIntron = new ChrBaseRegion(
                    candidate.chromosome(), up.RefEnd - upLoss + 1, down.Start + downLoss - 1);
            final int tier = classifyJunctionTier(candidateIntron);
            if(tier == SpliceMotif.TIER_NONE)
                continue;

            final int candidateMinAnchor = Math.min(up.TrailingM - upLoss, down.LeadingM - downLoss);
            if(tier > chosenTier || (tier == chosenTier && candidateMinAnchor > chosenMinAnchor))
            {
                chosenL = L;
                chosenIntron = candidateIntron;
                chosenTier = tier;
                chosenMinAnchor = candidateMinAnchor;
            }
        }
        return new SnapPick(chosenL, chosenIntron);
    }

    // Use the mate's previously-rescued intron as a hint. Primary-upstream pins intron start;
    // primary-downstream pins intron end.
    private SnapPick scanMateHint(
            final RescueCandidate candidate, final Side up, final Side down,
            final int upMatchedRead, final boolean primaryIsUpstream)
    {
        if(candidate.mateHintIntrons().isEmpty())
            return SnapPick.miss();

        for(ChrBaseRegion hint : candidate.mateHintIntrons())
        {
            if(!hint.Chromosome.equals(candidate.chromosome()))
                continue;
            final int upLoss;
            if(primaryIsUpstream)
                upLoss = up.RefEnd - hint.start() + 1;
            else
                upLoss = upMatchedRead - (down.LeadingS + (hint.end() - (down.Start - 1)));
            final int L = upMatchedRead - upLoss;
            if(L < down.LeadingS || L > upMatchedRead)
                continue;
            final int downLoss = L - down.LeadingS;
            if(up.TrailingM - upLoss < mConfig.MinAnchorOverhang)
                continue;
            if(down.LeadingM - downLoss < mConfig.MinAnchorOverhang)
                continue;
            if(up.TrailingM < upLoss || down.LeadingM < downLoss)
                continue;

            final int hintedIntronStart = primaryIsUpstream ? hint.start() : (up.RefEnd - upLoss + 1);
            final int hintedIntronEnd = primaryIsUpstream ? (down.Start + downLoss - 1) : hint.end();
            return new SnapPick(L, new ChrBaseRegion(candidate.chromosome(), hintedIntronStart, hintedIntronEnd));
        }
        return SnapPick.miss();
    }

    private static List<CigarElement> buildMergedCigar(
            final List<CigarElement> upCigar, final List<CigarElement> downCigar,
            final int upLoss, final int downLoss, final int intronLength)
    {
        final List<CigarElement> merged = new ArrayList<>(upCigar.size() + downCigar.size());
        for(int i = 0; i < upCigar.size() - 1; ++i) // upstream ops, excluding trailing S
        {
            if(i == upCigar.size() - 2 && upLoss > 0)
                merged.add(new CigarElement(upCigar.get(i).getLength() - upLoss, upCigar.get(i).getOperator()));
            else
                merged.add(upCigar.get(i));
        }
        merged.add(new CigarElement(intronLength, CigarOperator.N));
        for(int i = 1; i < downCigar.size(); ++i) // downstream ops, excluding leading S
        {
            if(i == 1 && downLoss > 0)
                merged.add(new CigarElement(downCigar.get(i).getLength() - downLoss, downCigar.get(i).getOperator()));
            else
                merged.add(downCigar.get(i));
        }
        return merged;
    }

    // Length of the M/=/X run adjacent to the trailing softclip (e.g. "5M 50S" -> 5).
    private static int trailingMatchedRun(final List<CigarElement> elements)
    {
        for(int i = elements.size() - 1; i >= 0; --i)
        {
            final CigarOperator op = elements.get(i).getOperator();
            if(op == CigarOperator.S)
                continue;
            if(op == CigarOperator.M || op == CigarOperator.EQ || op == CigarOperator.X)
                return elements.get(i).getLength();
            return 0;
        }
        return 0;
    }

    private static int leadingMatchedRun(final List<CigarElement> elements)
    {
        for(int i = 0; i < elements.size(); ++i)
        {
            final CigarOperator op = elements.get(i).getOperator();
            if(op == CigarOperator.S)
                continue;
            if(op == CigarOperator.M || op == CigarOperator.EQ || op == CigarOperator.X)
                return elements.get(i).getLength();
            return 0;
        }
        return 0;
    }

    // Folds tx-contig lift artifact: "<block>M nN yM zS" (trailing) or mirror (leading), where yM
    // is sub-threshold, into "<block>M (y+z)S". Drops the spurious N and conserves read-base count.
    static List<CigarElement> foldFabricatedTerminalMicroJunctions(
            final List<CigarElement> elements, final int minAnchor)
    {
        List<CigarElement> result = foldTrailingFabricatedMicroJunction(elements, minAnchor);
        result = foldLeadingFabricatedMicroJunction(result, minAnchor);
        return result;
    }

    private static List<CigarElement> foldTrailingFabricatedMicroJunction(
            final List<CigarElement> elements, final int minAnchor)
    {
        final int n = elements.size();
        if(n < 4)
            return elements;
        if(elements.get(n - 1).getOperator() != CigarOperator.S)
            return elements;
        final CigarElement anchor = elements.get(n - 2);
        if(anchor.getOperator() != CigarOperator.M || anchor.getLength() >= minAnchor)
            return elements;
        if(elements.get(n - 3).getOperator() != CigarOperator.N)
            return elements;
        if(!isMatchedOp(elements.get(n - 4).getOperator())) // real matched block required before the N
            return elements;

        final int foldedClip = anchor.getLength() + elements.get(n - 1).getLength();
        final List<CigarElement> out = new ArrayList<>(elements.subList(0, n - 3));
        out.add(new CigarElement(foldedClip, CigarOperator.S));
        return out;
    }

    private static List<CigarElement> foldLeadingFabricatedMicroJunction(
            final List<CigarElement> elements, final int minAnchor)
    {
        final int n = elements.size();
        if(n < 4)
            return elements;
        if(elements.get(0).getOperator() != CigarOperator.S)
            return elements;
        final CigarElement anchor = elements.get(1);
        if(anchor.getOperator() != CigarOperator.M || anchor.getLength() >= minAnchor)
            return elements;
        if(elements.get(2).getOperator() != CigarOperator.N)
            return elements;
        if(!isMatchedOp(elements.get(3).getOperator())) // real matched block required after the N
            return elements;

        final int foldedClip = elements.get(0).getLength() + anchor.getLength();
        final List<CigarElement> out = new ArrayList<>(n - 2);
        out.add(new CigarElement(foldedClip, CigarOperator.S));
        out.addAll(elements.subList(3, n));
        return out;
    }

    private static boolean isMatchedOp(final CigarOperator op)
    {
        return op == CigarOperator.M || op == CigarOperator.EQ || op == CigarOperator.X;
    }

    private static boolean opAdjacentToSoftClip(final List<CigarElement> elements, final boolean leadingSide)
    {
        if(elements.size() < 2)
            return false;
        if(leadingSide)
        {
            if(elements.get(0).getOperator() != CigarOperator.S)
                return false;
            final CigarOperator next = elements.get(1).getOperator();
            return next == CigarOperator.I || next == CigarOperator.D;
        }
        else
        {
            final int last = elements.size() - 1;
            if(elements.get(last).getOperator() != CigarOperator.S)
                return false;
            final CigarOperator prev = elements.get(last - 1).getOperator();
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
            return null;

        final int suppRefEnd = suppStart + CigarUtils.cigarAlignedLength(suppCigar) - 1;
        if(suppStart < primaryStart)
            return clampUpstreamSupp(suppCigar, suppStart, primaryStart);
        if(suppRefEnd > primaryRefEnd)
            return clampDownstreamSupp(suppCigar, suppStart, primaryRefEnd);
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
            final CigarElement e = suppCigar.get(i);
            final CigarOperator op = e.getOperator();
            final int len = e.getLength();

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
                refCursor += len;
            if(op.consumesReadBases())
                readCursor += len;
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
            final CigarElement e = suppCigar.get(i);
            final CigarOperator op = e.getOperator();
            final int len = e.getLength();

            if(op == CigarOperator.M || op == CigarOperator.EQ || op == CigarOperator.X)
            {
                final int mRefEnd = refCursor + len - 1;
                prevMInsidePrimary = mRefEnd <= primaryRefEnd;
            }

            if(op == CigarOperator.N && prevMInsidePrimary)
            {
                lastBadNIdx = i;
                refAfterLastBadN = refCursor + len; // just after this N op
                readAfterLastBadN = readCursor;
            }

            if(op.consumesReferenceBases())
                refCursor += len;
            if(op.consumesReadBases())
                readCursor += len;
        }

        if(lastBadNIdx == -1 || refAfterLastBadN < 0)
            return null;
        return cutHeadAt(suppCigar, refAfterLastBadN, lastBadNIdx, readAfterLastBadN);
    }

    private static ClampedSupp cutTailAt(
            final List<CigarElement> suppCigar, final int suppStart,
            final int cutIdx, final int readBeforeCut)
    {
        final List<CigarElement> trimmed = new ArrayList<>(cutIdx + 1);
        for(int i = 0; i < cutIdx; i++)
            trimmed.add(suppCigar.get(i));
        final int totalRead = CigarUtils.cigarBaseLength(suppCigar);
        final int trailingS = totalRead - readBeforeCut;
        if(trailingS > 0)
            trimmed.add(new CigarElement(trailingS, CigarOperator.S));
        return new ClampedSupp(suppStart, trimmed);
    }

    private static ClampedSupp cutHeadAt(
            final List<CigarElement> suppCigar, final int newStart,
            final int cutIdx, final int readBeforeCut)
    {
        final List<CigarElement> trimmed = new ArrayList<>(suppCigar.size() - cutIdx);
        if(readBeforeCut > 0)
            trimmed.add(new CigarElement(readBeforeCut, CigarOperator.S));
        for(int i = cutIdx + 1; i < suppCigar.size(); i++)
            trimmed.add(suppCigar.get(i));
        return new ClampedSupp(newStart, trimmed);
    }

    private static final class ClampedSupp
    {
        final int Start;
        final List<CigarElement> Cigar;

        ClampedSupp(final int start, final List<CigarElement> cigar)
        {
            Start = start;
            Cigar = cigar;
        }
    }

    private static final class Side
    {
        final int Start;
        final List<CigarElement> Cigar;
        final int LeadingS;
        final int TrailingS;
        final int LeadingM;
        final int TrailingM;
        final int RefEnd;

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
        final int L;
        final ChrBaseRegion Intron;

        SnapPick(final int l, final ChrBaseRegion intron)
        {
            L = l;
            Intron = intron;
        }

        static SnapPick miss()
        {
            return new SnapPick(-1, null);
        }
    }

    private static final class MergeOutcome
    {
        final RescueRejectReason Reject;
        final int MergedStart;
        final List<CigarElement> MergedCigar;
        final ChrBaseRegion IntroducedIntron;
        final RescueSupplementary MergedSupp;
        final boolean RightExtend;       // which terminal softclip of the primary this merge resolved

        private MergeOutcome(
                final RescueRejectReason reject,
                final int mergedStart, final List<CigarElement> mergedCigar,
                final ChrBaseRegion intron, final RescueSupplementary supp, final boolean rightExtend)
        {
            Reject = reject;
            MergedStart = mergedStart;
            MergedCigar = mergedCigar;
            IntroducedIntron = intron;
            MergedSupp = supp;
            RightExtend = rightExtend;
        }

        static MergeOutcome reject(final RescueRejectReason reason)
        {
            return new MergeOutcome(reason, -1, null, null, null, false);
        }

        static MergeOutcome success(
                final int start, final List<CigarElement> cigar,
                final ChrBaseRegion intron, final RescueSupplementary supp, final boolean rightExtend)
        {
            return new MergeOutcome(null, start, cigar, intron, supp, rightExtend);
        }

        boolean isSuccess()
        {
            return Reject == null;
        }
    }
}
