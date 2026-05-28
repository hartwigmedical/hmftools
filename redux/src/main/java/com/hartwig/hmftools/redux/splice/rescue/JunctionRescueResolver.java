package com.hartwig.hmftools.redux.splice.rescue;

import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_DELETION;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_INSERTION;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_MATCH;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_SEQ_MATCH;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_SEQ_MISMATCH;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_SKIPPED;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_SOFTCLIP;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

// Detects when a primary record's terminal softclip is actually a splice across an annotated
// intron, with the missing exon-side already emitted by bwa as a supplementary record. Merges the
// primary + supp into a single spliced primary with an N op.
//
// Per-mate operation: caller assembles a RescueCandidate (primary + supps list) and calls resolve().
// Result tells caller (a) whether the primary was rewritten and (b) which supp indices to drop.
// The resolver chains merges up to RescueConfig.MaxChainDepth so a read split into 3 alignments
// (one primary + two supps spanning 3 exons) can collapse into a single spliced primary.
public class JunctionRescueResolver
{
    private final AnnotatedJunctionIndex mAnnotatedIndex;
    private final RefSequenceSource mRefSource;
    private final RescueConfig mConfig;
    private final RescueStatistics mStatistics;

    public JunctionRescueResolver(final Set<ChrIntron> annotatedJunctions, final RescueConfig config)
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

        if(candidate.Supplementaries.isEmpty())
            return tryRefVerifyOnly(candidate);

        mStatistics.countCandidate();

        List<CigarShape.Element> primaryCigar = CigarShape.parse(candidate.PrimaryCigar);
        if(CigarShape.hasHardClip(primaryCigar))
        {
            mStatistics.countReject(RescueRejectReason.COMPLEX_CIGAR_SHAPE);
            return RescueResult.noMerge(RescueRejectReason.COMPLEX_CIGAR_SHAPE);
        }

        int primaryStart = candidate.PrimaryStart;
        final List<RescueSupplementary> remaining = new ArrayList<>(candidate.Supplementaries);
        final List<Integer> dropped = new ArrayList<>();
        final List<ChrIntron> introns = new ArrayList<>();
        RescueRejectReason lastReject = null;
        int chainDepth = 0;

        while(chainDepth < mConfig.MaxChainDepth && !remaining.isEmpty())
        {
            final boolean primaryHasLeadingS = !primaryCigar.isEmpty()
                    && primaryCigar.get(0).Op == OP_SOFTCLIP;
            final boolean primaryHasTrailingS = !primaryCigar.isEmpty()
                    && primaryCigar.get(primaryCigar.size() - 1).Op == OP_SOFTCLIP;
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
            dropped.add(attempt.MergedSupp.Index);
            introns.add(attempt.IntroducedIntron);
            remaining.remove(attempt.MergedSupp);
            ++chainDepth;
        }

        if(chainDepth == 0)
        {
            // shouldn't reach here — initial-iteration failures returned above. Belt and braces.
            return RescueResult.noMerge(lastReject != null ? lastReject : RescueRejectReason.NO_MATCHING_SUPP);
        }

        mStatistics.countMergedChain(chainDepth);
        return new RescueResult(
                true, CigarShape.format(primaryCigar), primaryStart,
                dropped, introns, chainDepth, null);
    }

    // Ref-verify path: primary has a terminal softclip but no supplementary. Look up annotated
    // junctions whose start/end abuts the softclip boundary, fetch reference bases at the
    // candidate next-exon position, and accept if the softclipped read bases match within tolerance.
    private RescueResult tryRefVerifyOnly(final RescueCandidate candidate)
    {
        if(mRefSource == null || candidate.ReadBases == null)
            return RescueResult.noMerge(RescueRejectReason.NO_MATCHING_SUPP);

        mStatistics.countCandidate();

        final List<CigarShape.Element> primaryCigar = CigarShape.parse(candidate.PrimaryCigar);
        if(CigarShape.hasHardClip(primaryCigar))
        {
            mStatistics.countReject(RescueRejectReason.COMPLEX_CIGAR_SHAPE);
            return RescueResult.noMerge(RescueRejectReason.COMPLEX_CIGAR_SHAPE);
        }
        if(primaryCigar.isEmpty())
            return RescueResult.noMerge(RescueRejectReason.NO_TERMINAL_SOFTCLIP);

        final boolean trailingS = primaryCigar.get(primaryCigar.size() - 1).Op == OP_SOFTCLIP;
        final boolean leadingS = primaryCigar.get(0).Op == OP_SOFTCLIP;
        if(!trailingS && !leadingS)
            return RescueResult.noMerge(RescueRejectReason.NO_TERMINAL_SOFTCLIP);
        if(trailingS && leadingS)
            return RescueResult.noMerge(RescueRejectReason.COMPLEX_CIGAR_SHAPE);

        final boolean rightExtend = trailingS;
        if(opAdjacentToSoftClip(primaryCigar, !rightExtend))
            return RescueResult.noMerge(RescueRejectReason.COMPLEX_CIGAR_SHAPE);

        final int softclipLen = rightExtend
                ? CigarShape.trailingSoftClip(primaryCigar)
                : CigarShape.leadingSoftClip(primaryCigar);
        if(softclipLen < mConfig.MinAnchorOverhang)
            return RescueResult.noMerge(RescueRejectReason.SHORT_ANCHOR);
        final int primaryAnchor = rightExtend
                ? trailingMatchedRun(primaryCigar)
                : leadingMatchedRun(primaryCigar);
        if(primaryAnchor < mConfig.MinAnchorOverhang)
            return RescueResult.noMerge(RescueRejectReason.SHORT_ANCHOR);

        final List<ChrIntron> candidates;
        if(rightExtend)
        {
            final int primaryRefEnd = candidate.PrimaryStart + CigarShape.referenceSpan(primaryCigar) - 1;
            candidates = mAnnotatedIndex.introByStart(candidate.Chromosome, primaryRefEnd + 1);
        }
        else
        {
            candidates = mAnnotatedIndex.introByEnd(candidate.Chromosome, candidate.PrimaryStart - 1);
        }
        return verifyAgainstCandidates(candidate, primaryCigar, candidates, softclipLen, rightExtend);
    }

    private RescueResult verifyAgainstCandidates(
            final RescueCandidate candidate, final List<CigarShape.Element> primaryCigar,
            final List<ChrIntron> candidates, final int softclipLen, final boolean rightExtend)
    {
        if(candidates.isEmpty())
        {
            mStatistics.countReject(RescueRejectReason.REF_VERIFY_NO_CANDIDATE_EXON);
            return RescueResult.noMerge(RescueRejectReason.REF_VERIFY_NO_CANDIDATE_EXON);
        }

        final byte[] readBases = candidate.ReadBases;
        final byte[] softclipBases = new byte[softclipLen];
        if(rightExtend)
            System.arraycopy(readBases, readBases.length - softclipLen, softclipBases, 0, softclipLen);
        else
            System.arraycopy(readBases, 0, softclipBases, 0, softclipLen);

        final int maxMismatches = softclipLen / 10;        // ~90% identity floor

        ChrIntron chosen = null;
        int chosenMismatches = Integer.MAX_VALUE;
        int ambiguous = 0;

        for(ChrIntron candidateIntron : candidates)
        {
            final int intronLength = candidateIntron.IntronEnd - candidateIntron.IntronStart + 1;
            if(intronLength < mConfig.MinIntronLength)
                continue;
            if(intronLength > mConfig.MaxIntronLength)
                continue;

            final int refStart;
            final int refEnd;
            if(rightExtend)
            {
                refStart = candidateIntron.IntronEnd + 1;
                refEnd = candidateIntron.IntronEnd + softclipLen;
            }
            else
            {
                refStart = candidateIntron.IntronStart - softclipLen;
                refEnd = candidateIntron.IntronStart - 1;
            }

            final byte[] refBases = mRefSource.getBases(candidate.Chromosome, refStart, refEnd);
            if(refBases == null || refBases.length != softclipLen)
                continue;

            int mismatches = 0;
            for(int i = 0; i < softclipLen && mismatches <= maxMismatches; ++i)
            {
                if(!basesEqualIgnoreCase(softclipBases[i], refBases[i]))
                    ++mismatches;
            }
            if(mismatches > maxMismatches)
                continue;

            if(chosen == null || mismatches < chosenMismatches)
            {
                chosen = candidateIntron;
                chosenMismatches = mismatches;
                ambiguous = 0;
            }
            else if(mismatches == chosenMismatches)
            {
                ++ambiguous;
            }
        }

        if(chosen == null)
        {
            mStatistics.countReject(RescueRejectReason.REF_VERIFY_MISMATCH_TOO_HIGH);
            return RescueResult.noMerge(RescueRejectReason.REF_VERIFY_MISMATCH_TOO_HIGH);
        }
        if(ambiguous > 0)
        {
            mStatistics.countReject(RescueRejectReason.REF_VERIFY_AMBIGUOUS);
            return RescueResult.noMerge(RescueRejectReason.REF_VERIFY_AMBIGUOUS);
        }

        final int intronLength = chosen.IntronEnd - chosen.IntronStart + 1;
        final List<CigarShape.Element> merged = new ArrayList<>(primaryCigar.size() + 2);
        if(rightExtend)
        {
            for(int i = 0; i < primaryCigar.size() - 1; ++i)
                merged.add(primaryCigar.get(i));
            merged.add(new CigarShape.Element(intronLength, OP_SKIPPED));
            merged.add(new CigarShape.Element(softclipLen, OP_MATCH));
        }
        else
        {
            merged.add(new CigarShape.Element(softclipLen, OP_MATCH));
            merged.add(new CigarShape.Element(intronLength, OP_SKIPPED));
            for(int i = 1; i < primaryCigar.size(); ++i)
                merged.add(primaryCigar.get(i));
        }

        final int mergedStart = rightExtend ? candidate.PrimaryStart : (chosen.IntronStart - softclipLen);
        mStatistics.countMergedChain(1);
        return new RescueResult(
                true, CigarShape.format(merged), mergedStart,
                Collections.emptyList(), Collections.singletonList(chosen),
                1, null);
    }

    // Annotation > motif. Returns TIER_NONE when ref-source is unavailable.
    private int classifyJunctionTier(final ChrIntron candidateIntron)
    {
        if(mAnnotatedIndex.contains(candidateIntron))
            return SpliceMotif.TIER_ANNOTATED;
        if(mRefSource == null)
            return SpliceMotif.TIER_NONE;
        final byte[] donor = mRefSource.getBases(
                candidateIntron.Chromosome, candidateIntron.IntronStart, candidateIntron.IntronStart + 1);
        final byte[] acceptor = mRefSource.getBases(
                candidateIntron.Chromosome, candidateIntron.IntronEnd - 1, candidateIntron.IntronEnd);
        return SpliceMotif.classify(donor, acceptor);
    }

    private static boolean basesEqualIgnoreCase(final byte a, final byte b)
    {
        if(a == b)
            return true;
        // mask the 0x20 bit (ASCII case bit) and re-compare; only valid for letters but cheap.
        return (a & ~0x20) == (b & ~0x20);
    }

    // Picks the highest-confidence supplementary that passes all gates. Returns the merged-primary
    // proposal or a reject reason describing which gate the best candidate(s) failed against.
    private MergeOutcome pickBestSupplementary(
            final RescueCandidate candidate, final int primaryStart,
            final List<CigarShape.Element> primaryCigar, final List<RescueSupplementary> supps)
    {
        MergeOutcome best = null;
        RescueRejectReason lastReject = null;

        for(RescueSupplementary supp : supps)
        {
            final MergeOutcome attempt = tryMerge(candidate, primaryStart, primaryCigar, supp);
            if(!attempt.isSuccess())
            {
                mStatistics.countReject(attempt.Reject);
                lastReject = attempt.Reject;
                continue;
            }

            if(best == null)
            {
                best = attempt;
            }
            else if(best.Ambiguous)
            {
                // already ambiguous; a third equally-good candidate doesn't recover. Keep iterating
                // only to drain reject stats on remaining failures via the !isSuccess branch above.
            }
            else if(isBetterCandidate(attempt, best))
            {
                best = attempt;
            }
            else if(isTied(attempt, best))
            {
                best = MergeOutcome.ambiguous();
            }
        }

        if(best == null)
            return MergeOutcome.reject(lastReject != null ? lastReject : RescueRejectReason.NO_MATCHING_SUPP);

        if(best.Ambiguous)
        {
            mStatistics.countReject(RescueRejectReason.AMBIGUOUS_SUPP_CHOICE);
            return MergeOutcome.reject(RescueRejectReason.AMBIGUOUS_SUPP_CHOICE);
        }

        return best;
    }

    // Candidate ranking. Axes in order: higher MAPQ wins, then smaller intron wins. A tie on both
    // becomes AMBIGUOUS_SUPP_CHOICE (caller refuses to merge).
    private static boolean isBetterCandidate(final MergeOutcome a, final MergeOutcome b)
    {
        if(a.MergedSupp.Mapq != b.MergedSupp.Mapq)
            return a.MergedSupp.Mapq > b.MergedSupp.Mapq;
        return intronLength(a.IntroducedIntron) < intronLength(b.IntroducedIntron);
    }

    private static boolean isTied(final MergeOutcome a, final MergeOutcome b)
    {
        return a.MergedSupp.Mapq == b.MergedSupp.Mapq
                && intronLength(a.IntroducedIntron) == intronLength(b.IntroducedIntron);
    }

    private static int intronLength(final ChrIntron intron)
    {
        return intron.IntronEnd - intron.IntronStart + 1;
    }

    private MergeOutcome tryMerge(
            final RescueCandidate candidate, final int primaryStart,
            final List<CigarShape.Element> primaryCigar, final RescueSupplementary supp)
    {
        if(!candidate.Chromosome.equals(supp.Chromosome))
            return MergeOutcome.reject(RescueRejectReason.DIFFERENT_CHROMOSOME);

        if(candidate.ForwardStrand != supp.ForwardStrand)
            return MergeOutcome.reject(RescueRejectReason.OPPOSITE_STRAND);

        final List<CigarShape.Element> suppCigar = CigarShape.parse(supp.Cigar);
        if(CigarShape.hasHardClip(suppCigar))
            return MergeOutcome.reject(RescueRejectReason.COMPLEX_CIGAR_SHAPE);

        // Decide which record is upstream and which is downstream by softclip shapes (and, when a
        // middle-anchored primary has S on both sides, by where the supp sits genomically).
        final int primaryLeadS = CigarShape.leadingSoftClip(primaryCigar);
        final int primaryTrailS = CigarShape.trailingSoftClip(primaryCigar);
        final int suppLeadS = CigarShape.leadingSoftClip(suppCigar);
        final int suppTrailS = CigarShape.trailingSoftClip(suppCigar);

        boolean rightExtend = primaryTrailS > 0 && suppLeadS > 0;
        boolean leftExtend = primaryLeadS > 0 && suppTrailS > 0;

        if(!rightExtend && !leftExtend)
            return MergeOutcome.reject(RescueRejectReason.NO_MATCHING_SUPP);

        if(rightExtend && leftExtend)
        {
            final int primaryRefEnd = primaryStart + CigarShape.referenceSpan(primaryCigar) - 1;
            final int suppRefEnd = supp.Start + CigarShape.referenceSpan(suppCigar) - 1;
            if(supp.Start > primaryRefEnd && suppRefEnd >= primaryStart)
                leftExtend = false;
            else if(suppRefEnd < primaryStart && supp.Start <= primaryRefEnd)
                rightExtend = false;
            else
                return MergeOutcome.reject(RescueRejectReason.COMPLEX_CIGAR_SHAPE);
        }

        final Side primarySide = Side.of(primaryStart, primaryCigar);
        final Side suppSide = Side.of(supp.Start, suppCigar);
        final boolean primaryIsUpstream = rightExtend;
        final Side up = primaryIsUpstream ? primarySide : suppSide;
        final Side down = primaryIsUpstream ? suppSide : primarySide;

        return mergeJunction(candidate, up, down, primaryIsUpstream, supp);
    }

    // The unified merge: validates the up/down anchor pair, scores candidate snap points by
    // junction tier and min anchor, falls back to mate-hint, then to trust-primary. The body is
    // direction-agnostic — caller has already mapped primary/supp into upstream/downstream slots.
    // supp is threaded through purely as the result payload (so callers know which supp was merged).
    private MergeOutcome mergeJunction(
            final RescueCandidate candidate, final Side up, final Side down,
            final boolean primaryIsUpstream, final RescueSupplementary supp)
    {
        if(CigarShape.readLength(up.Cigar) != candidate.ReadLength
                || CigarShape.readLength(down.Cigar) != candidate.ReadLength)
            return MergeOutcome.reject(RescueRejectReason.COMPLEX_CIGAR_SHAPE);

        if(opAdjacentToSoftClip(up.Cigar, false) || opAdjacentToSoftClip(down.Cigar, true))
            return MergeOutcome.reject(RescueRejectReason.COMPLEX_CIGAR_SHAPE);

        final int upMatchedRead = candidate.ReadLength - up.TrailingS;
        final int overlap = upMatchedRead - down.LeadingS;

        if(overlap < 0)
            return MergeOutcome.reject(RescueRejectReason.READ_COVERAGE_GAP);
        if(overlap > mConfig.SoftclipTolerance)
            return MergeOutcome.reject(RescueRejectReason.READ_COVERAGE_OVERLAP);
        if(down.Start <= up.RefEnd)
            return MergeOutcome.reject(RescueRejectReason.READ_COVERAGE_OVERLAP);

        // Intron length is invariant under the snap point L: any redistribution of overlap bases
        // moves both intron endpoints by the same amount. Check once, up front.
        final int intronLength = (down.Start - 1 - up.RefEnd) + overlap;
        if(intronLength < mConfig.MinIntronLength)
            return MergeOutcome.reject(RescueRejectReason.INTRON_TOO_SHORT);
        if(intronLength > mConfig.MaxIntronLength)
            return MergeOutcome.reject(RescueRejectReason.INTRON_TOO_LONG);

        if(up.TrailingM < mConfig.MinAnchorOverhang || down.LeadingM < mConfig.MinAnchorOverhang)
            return MergeOutcome.reject(RescueRejectReason.SHORT_ANCHOR);

        // Pick the highest-tier snap point in L ∈ [down.LeadingS, upMatchedRead]; within a tier,
        // largest min(upAnchor, downAnchor) wins; ties go to "first encountered". Iterate so the
        // primary-fully-matched value (upLoss=0 when primary is upstream, downLoss=0 when primary
        // is downstream) is visited first — that preserves trust-primary as the natural tie-break.
        final SnapPick pick = scanSnapPoints(candidate, up, down, upMatchedRead, primaryIsUpstream);
        int chosenL = pick.L;
        ChrIntron chosenIntron = pick.Intron;

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
            // Trust-primary fallback: primary contributes everything it matched. If primary is
            // upstream, upLoss=0 → L=upMatchedRead. If primary is downstream, downLoss=0 → L=down.LeadingS.
            chosenL = primaryIsUpstream ? upMatchedRead : down.LeadingS;
            final int upLossFb = upMatchedRead - chosenL;
            final int downLossFb = chosenL - down.LeadingS;
            if(up.TrailingM - upLossFb < mConfig.MinAnchorOverhang || up.TrailingM < upLossFb
                    || down.LeadingM - downLossFb < mConfig.MinAnchorOverhang || down.LeadingM < downLossFb)
                return MergeOutcome.reject(RescueRejectReason.SHORT_ANCHOR);
            chosenIntron = new ChrIntron(candidate.Chromosome, up.RefEnd + 1, down.Start - 1);
        }

        final int upLoss = upMatchedRead - chosenL;
        final int downLoss = chosenL - down.LeadingS;
        final List<CigarShape.Element> merged = buildMergedCigar(up.Cigar, down.Cigar, upLoss, downLoss, intronLength);
        return MergeOutcome.success(up.Start, merged, chosenIntron, supp);
    }

    private SnapPick scanSnapPoints(
            final RescueCandidate candidate, final Side up, final Side down,
            final int upMatchedRead, final boolean primaryIsUpstream)
    {
        int chosenL = -1;
        ChrIntron chosenIntron = null;
        int chosenTier = SpliceMotif.TIER_NONE;
        int chosenMinAnchor = -1;

        // Iterate so trust-primary is first encountered (wins ties).
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

            final ChrIntron candidateIntron = new ChrIntron(
                    candidate.Chromosome, up.RefEnd - upLoss + 1, down.Start + downLoss - 1);
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

    // Mate-hint fallback: partner mate's previously-rescued intron hints the junction position.
    // The hint position that's load-bearing depends on which side is the primary: if primary is
    // upstream, the hint pins the intron's start (up.RefEnd-side); if downstream, the intron's
    // end (down.Start-side).
    private SnapPick scanMateHint(
            final RescueCandidate candidate, final Side up, final Side down,
            final int upMatchedRead, final boolean primaryIsUpstream)
    {
        if(candidate.MateHintIntrons.isEmpty())
            return SnapPick.miss();

        for(ChrIntron hint : candidate.MateHintIntrons)
        {
            if(!hint.Chromosome.equals(candidate.Chromosome))
                continue;
            final int upLoss;
            if(primaryIsUpstream)
                upLoss = up.RefEnd - hint.IntronStart + 1;
            else
                upLoss = upMatchedRead - (down.LeadingS + (hint.IntronEnd - (down.Start - 1)));
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

            // Hint pins the primary's side of the intron; the other end is derived from the snap.
            final int hintedIntronStart = primaryIsUpstream ? hint.IntronStart : (up.RefEnd - upLoss + 1);
            final int hintedIntronEnd = primaryIsUpstream ? (down.Start + downLoss - 1) : hint.IntronEnd;
            return new SnapPick(L, new ChrIntron(candidate.Chromosome, hintedIntronStart, hintedIntronEnd));
        }
        return SnapPick.miss();
    }

    private static List<CigarShape.Element> buildMergedCigar(
            final List<CigarShape.Element> upCigar, final List<CigarShape.Element> downCigar,
            final int upLoss, final int downLoss, final int intronLength)
    {
        final List<CigarShape.Element> merged = new ArrayList<>(upCigar.size() + downCigar.size());
        // upstream ops, excluding trailing S; shrink last M by upLoss
        for(int i = 0; i < upCigar.size() - 1; ++i)
        {
            if(i == upCigar.size() - 2 && upLoss > 0)
                merged.add(new CigarShape.Element(upCigar.get(i).Length - upLoss, upCigar.get(i).Op));
            else
                merged.add(upCigar.get(i));
        }
        merged.add(new CigarShape.Element(intronLength, OP_SKIPPED));
        // downstream ops, excluding leading S; shrink first M by downLoss
        for(int i = 1; i < downCigar.size(); ++i)
        {
            if(i == 1 && downLoss > 0)
                merged.add(new CigarShape.Element(downCigar.get(i).Length - downLoss, downCigar.get(i).Op));
            else
                merged.add(downCigar.get(i));
        }
        return merged;
    }

    // returns the length of the matched-run (M/=/X) adjacent to the cigar's trailing softclip. If
    // the cigar ends in `... 5M 50S`, this returns 5.
    private static int trailingMatchedRun(final List<CigarShape.Element> elements)
    {
        for(int i = elements.size() - 1; i >= 0; --i)
        {
            final char op = elements.get(i).Op;
            if(op == OP_SOFTCLIP)
                continue;
            if(op == OP_MATCH || op == OP_SEQ_MATCH || op == OP_SEQ_MISMATCH)
                return elements.get(i).Length;
            return 0;
        }
        return 0;
    }

    private static int leadingMatchedRun(final List<CigarShape.Element> elements)
    {
        for(int i = 0; i < elements.size(); ++i)
        {
            final char op = elements.get(i).Op;
            if(op == OP_SOFTCLIP)
                continue;
            if(op == OP_MATCH || op == OP_SEQ_MATCH || op == OP_SEQ_MISMATCH)
                return elements.get(i).Length;
            return 0;
        }
        return 0;
    }

    // checks whether an I or D op sits immediately next to the leading or trailing softclip. Such
    // ops complicate the merge because the matched anchor isn't a clean M block; v1 rejects them.
    private static boolean opAdjacentToSoftClip(final List<CigarShape.Element> elements, final boolean leadingSide)
    {
        if(elements.size() < 2)
            return false;
        if(leadingSide)
        {
            if(elements.get(0).Op != OP_SOFTCLIP)
                return false;
            final char next = elements.get(1).Op;
            return next == OP_INSERTION || next == OP_DELETION;
        }
        else
        {
            final int last = elements.size() - 1;
            if(elements.get(last).Op != OP_SOFTCLIP)
                return false;
            final char prev = elements.get(last - 1).Op;
            return prev == OP_INSERTION || prev == OP_DELETION;
        }
    }

    // Per-record geometry computed once at the entry to tryMerge — start, parsed cigar, and the
    // S/M lengths on each side plus the reference end coordinate.
    private static final class Side
    {
        final int Start;
        final List<CigarShape.Element> Cigar;
        final int LeadingS;
        final int TrailingS;
        final int LeadingM;
        final int TrailingM;
        final int RefEnd;

        private Side(
                final int start, final List<CigarShape.Element> cigar,
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

        static Side of(final int start, final List<CigarShape.Element> cigar)
        {
            return new Side(start, cigar,
                    CigarShape.leadingSoftClip(cigar), CigarShape.trailingSoftClip(cigar),
                    leadingMatchedRun(cigar), trailingMatchedRun(cigar),
                    start + CigarShape.referenceSpan(cigar) - 1);
        }
    }

    private static final class SnapPick
    {
        final int L;
        final ChrIntron Intron;

        SnapPick(final int l, final ChrIntron intron)
        {
            L = l;
            Intron = intron;
        }

        static SnapPick miss()
        {
            return new SnapPick(-1, null);
        }
    }

    // Result of attempting to merge one supp into the primary. Three shapes: success (carries the
    // merged geometry plus the supp identity), reject (carries a reason), or ambiguous (two equally
    // good candidates — treated as a special reject by the caller).
    private static final class MergeOutcome
    {
        final RescueRejectReason Reject;
        final boolean Ambiguous;
        final int MergedStart;
        final List<CigarShape.Element> MergedCigar;
        final ChrIntron IntroducedIntron;
        final RescueSupplementary MergedSupp;

        private MergeOutcome(
                final RescueRejectReason reject, final boolean ambiguous,
                final int mergedStart, final List<CigarShape.Element> mergedCigar,
                final ChrIntron intron, final RescueSupplementary supp)
        {
            Reject = reject;
            Ambiguous = ambiguous;
            MergedStart = mergedStart;
            MergedCigar = mergedCigar;
            IntroducedIntron = intron;
            MergedSupp = supp;
        }

        static MergeOutcome reject(final RescueRejectReason reason)
        {
            return new MergeOutcome(reason, false, -1, null, null, null);
        }

        static MergeOutcome ambiguous()
        {
            return new MergeOutcome(RescueRejectReason.AMBIGUOUS_SUPP_CHOICE, true, -1, null, null, null);
        }

        static MergeOutcome success(
                final int start, final List<CigarShape.Element> cigar,
                final ChrIntron intron, final RescueSupplementary supp)
        {
            return new MergeOutcome(null, false, start, cigar, intron, supp);
        }

        boolean isSuccess()
        {
            return Reject == null;
        }
    }
}
