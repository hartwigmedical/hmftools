package com.hartwig.hmftools.redux.splice.rescue;

import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_DELETION;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_INSERTION;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_MATCH;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_SEQ_MATCH;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_SEQ_MISMATCH;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_SKIPPED;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_SOFTCLIP;

import java.util.ArrayList;
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

        // Empty supplementary list — only the ref-verify path can rescue. Fall through to it.
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

            final MergeAttempt attempt = pickBestSupplementary(
                    candidate, primaryStart, primaryCigar, remaining);

            if(!attempt.success())
            {
                if(chainDepth == 0)
                    return RescueResult.noMerge(attempt.RejectReason);
                lastReject = attempt.RejectReason;
                break;
            }

            primaryStart = attempt.Inner.MergedStart;
            primaryCigar = attempt.Inner.MergedCigar;
            dropped.add(attempt.MergedSupp.Index);
            introns.add(attempt.Inner.IntroducedIntron);
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
    // junctions that start (right-extend) or end (left-extend) where the softclip would extend,
    // fetch reference bases at the candidate next-exon position, and accept if the softclipped
    // read bases match within tolerance. Drops zero supps (there were none), introduces one
    // intron, no chain.
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

        if(trailingS && !leadingS)
            return tryRefVerifyRightExtend(candidate, primaryCigar);
        if(leadingS && !trailingS)
            return tryRefVerifyLeftExtend(candidate, primaryCigar);
        return RescueResult.noMerge(RescueRejectReason.COMPLEX_CIGAR_SHAPE);
    }

    private RescueResult tryRefVerifyRightExtend(
            final RescueCandidate candidate, final List<CigarShape.Element> primaryCigar)
    {
        if(opAdjacentToSoftClip(primaryCigar, false))
            return RescueResult.noMerge(RescueRejectReason.COMPLEX_CIGAR_SHAPE);

        final int trailS = CigarShape.trailingSoftClip(primaryCigar);
        if(trailS < mConfig.MinAnchorOverhang)
            return RescueResult.noMerge(RescueRejectReason.SHORT_ANCHOR);
        final int primaryAnchor = trailingMatchedRun(primaryCigar);
        if(primaryAnchor < mConfig.MinAnchorOverhang)
            return RescueResult.noMerge(RescueRejectReason.SHORT_ANCHOR);

        final int primaryRefEnd = candidate.PrimaryStart
                + CigarShape.referenceSpan(primaryCigar) - 1;
        final int intronStart = primaryRefEnd + 1;
        final List<ChrIntron> candidates = mAnnotatedIndex.introByStart(candidate.Chromosome, intronStart);
        return verifyAgainstCandidates(candidate, primaryCigar, candidates, trailS, true);
    }

    private RescueResult tryRefVerifyLeftExtend(
            final RescueCandidate candidate, final List<CigarShape.Element> primaryCigar)
    {
        if(opAdjacentToSoftClip(primaryCigar, true))
            return RescueResult.noMerge(RescueRejectReason.COMPLEX_CIGAR_SHAPE);

        final int leadS = CigarShape.leadingSoftClip(primaryCigar);
        if(leadS < mConfig.MinAnchorOverhang)
            return RescueResult.noMerge(RescueRejectReason.SHORT_ANCHOR);
        final int primaryAnchor = leadingMatchedRun(primaryCigar);
        if(primaryAnchor < mConfig.MinAnchorOverhang)
            return RescueResult.noMerge(RescueRejectReason.SHORT_ANCHOR);

        final int intronEnd = candidate.PrimaryStart - 1;
        final List<ChrIntron> candidates = mAnnotatedIndex.introByEnd(candidate.Chromosome, intronEnd);
        return verifyAgainstCandidates(candidate, primaryCigar, candidates, leadS, false);
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

        // The softclipped read bases: right-extend = trailing softclip's bases; left-extend = leading.
        final byte[] readBases = candidate.ReadBases;
        final byte[] softclipBases;
        if(rightExtend)
        {
            softclipBases = new byte[softclipLen];
            System.arraycopy(readBases, readBases.length - softclipLen, softclipBases, 0, softclipLen);
        }
        else
        {
            softclipBases = new byte[softclipLen];
            System.arraycopy(readBases, 0, softclipBases, 0, softclipLen);
        }

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
                // downstream exon's first softclipLen bases live at [intronEnd+1, intronEnd+softclipLen]
                refStart = candidateIntron.IntronEnd + 1;
                refEnd = candidateIntron.IntronEnd + softclipLen;
            }
            else
            {
                // upstream exon's last softclipLen bases end at intronStart-1
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

        // Rewrite cigar: drop terminal S, append N(intronLength) + softclipLen M.
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
                java.util.Collections.emptyList(), java.util.Collections.singletonList(chosen),
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
    private MergeAttempt pickBestSupplementary(
            final RescueCandidate candidate, final int primaryStart,
            final List<CigarShape.Element> primaryCigar, final List<RescueSupplementary> supps)
    {
        MergeAttempt bestSuccess = null;
        RescueAttemptResult lastAnyGateReject = null;

        for(RescueSupplementary supp : supps)
        {
            final RescueAttemptResult attempt = tryMerge(candidate, primaryStart, primaryCigar, supp);
            if(attempt.Reject != null)
            {
                mStatistics.countReject(attempt.Reject);
                lastAnyGateReject = attempt;
                continue;
            }

            if(bestSuccess == null)
            {
                bestSuccess = new MergeAttempt(attempt, supp);
                continue;
            }

            final int cmp = compareSuccess(attempt, supp, bestSuccess);
            if(cmp < 0)
            {
                bestSuccess = new MergeAttempt(attempt, supp);
            }
            else if(cmp == 0)
            {
                bestSuccess = MergeAttempt.ambiguous();
            }
        }

        if(bestSuccess == null)
        {
            final RescueRejectReason reason = lastAnyGateReject != null
                    ? lastAnyGateReject.Reject
                    : RescueRejectReason.NO_MATCHING_SUPP;
            return MergeAttempt.fail(reason);
        }

        if(bestSuccess.Ambiguous)
        {
            mStatistics.countReject(RescueRejectReason.AMBIGUOUS_SUPP_CHOICE);
            return MergeAttempt.fail(RescueRejectReason.AMBIGUOUS_SUPP_CHOICE);
        }

        return bestSuccess;
    }

    // Lower compareSuccess result == better candidate. Comparison axes (in order):
    //   1. higher MAPQ wins
    //   2. smaller intron length wins
    // Equal on both → ambiguous, caller refuses to merge.
    private int compareSuccess(final RescueAttemptResult a, final RescueSupplementary aSupp,
            final MergeAttempt current)
    {
        if(aSupp.Mapq != current.MergedSupp.Mapq)
            return current.MergedSupp.Mapq - aSupp.Mapq;       // higher MAPQ = "smaller" (better)
        final int aIntron = a.IntroducedIntron.IntronEnd - a.IntroducedIntron.IntronStart + 1;
        final ChrIntron curIntron = current.Inner.IntroducedIntron;
        return aIntron - (curIntron.IntronEnd - curIntron.IntronStart + 1);
    }

    private RescueAttemptResult tryMerge(
            final RescueCandidate candidate, final int primaryStart,
            final List<CigarShape.Element> primaryCigar, final RescueSupplementary supp)
    {
        if(!candidate.Chromosome.equals(supp.Chromosome))
            return RescueAttemptResult.reject(RescueRejectReason.DIFFERENT_CHROMOSOME);

        if(candidate.ForwardStrand != supp.ForwardStrand)
            return RescueAttemptResult.reject(RescueRejectReason.OPPOSITE_STRAND);

        final List<CigarShape.Element> suppCigar = CigarShape.parse(supp.Cigar);
        if(CigarShape.hasHardClip(suppCigar))
            return RescueAttemptResult.reject(RescueRejectReason.COMPLEX_CIGAR_SHAPE);

        // bail early if the cigar shapes don't make a clean right-extend or left-extend pair.
        final int primaryLeadS = CigarShape.leadingSoftClip(primaryCigar);
        final int primaryTrailS = CigarShape.trailingSoftClip(primaryCigar);
        final int suppLeadS = CigarShape.leadingSoftClip(suppCigar);
        final int suppTrailS = CigarShape.trailingSoftClip(suppCigar);

        // A middle-anchored primary in a 3-exon read has both leading AND trailing softclip and
        // chain-merges twice (right supp first, then left supp). Disambiguate by where the supp
        // sits genomically when shape alone matches both directions.
        boolean rightExtend = primaryTrailS > 0 && suppLeadS > 0;
        boolean leftExtend = primaryLeadS > 0 && suppTrailS > 0;

        if(!rightExtend && !leftExtend)
            return RescueAttemptResult.reject(RescueRejectReason.NO_MATCHING_SUPP);

        if(rightExtend && leftExtend)
        {
            final int primaryRefEnd = primaryStart + CigarShape.referenceSpan(primaryCigar) - 1;
            final int suppRefEnd = supp.Start + CigarShape.referenceSpan(suppCigar) - 1;
            if(supp.Start > primaryRefEnd && suppRefEnd >= primaryStart)
                leftExtend = false;
            else if(suppRefEnd < primaryStart && supp.Start <= primaryRefEnd)
                rightExtend = false;
            else
                return RescueAttemptResult.reject(RescueRejectReason.COMPLEX_CIGAR_SHAPE);
        }

        return rightExtend
                ? tryRightExtend(candidate, primaryStart, primaryCigar, supp, suppCigar)
                : tryLeftExtend(candidate, primaryStart, primaryCigar, supp, suppCigar);
    }

    private RescueAttemptResult tryRightExtend(
            final RescueCandidate candidate, final int primaryStart,
            final List<CigarShape.Element> primaryCigar, final RescueSupplementary supp,
            final List<CigarShape.Element> suppCigar)
    {
        final int primaryReadCovered = CigarShape.readLength(primaryCigar);
        final int suppReadCovered = CigarShape.readLength(suppCigar);

        if(primaryReadCovered != candidate.ReadLength || suppReadCovered != candidate.ReadLength)
            return RescueAttemptResult.reject(RescueRejectReason.COMPLEX_CIGAR_SHAPE);

        if(opAdjacentToSoftClip(primaryCigar, false))
            return RescueAttemptResult.reject(RescueRejectReason.COMPLEX_CIGAR_SHAPE);
        if(opAdjacentToSoftClip(suppCigar, true))
            return RescueAttemptResult.reject(RescueRejectReason.COMPLEX_CIGAR_SHAPE);

        final int primaryMatched = candidate.ReadLength - CigarShape.trailingSoftClip(primaryCigar);
        final int suppLeadingS = CigarShape.leadingSoftClip(suppCigar);
        final int overlap = primaryMatched - suppLeadingS;

        if(overlap < 0)
            return RescueAttemptResult.reject(RescueRejectReason.READ_COVERAGE_GAP);
        if(overlap > mConfig.SoftclipTolerance)
            return RescueAttemptResult.reject(RescueRejectReason.READ_COVERAGE_OVERLAP);

        final int primaryRefEnd = primaryStart + CigarShape.referenceSpan(primaryCigar) - 1;
        if(supp.Start <= primaryRefEnd)
            return RescueAttemptResult.reject(RescueRejectReason.READ_COVERAGE_OVERLAP);

        // Intron length is invariant under the snap point L (any redistribution of overlap bases
        // moves both intron endpoints by the same amount). Check once, up front.
        final int intronLength = (supp.Start - 1 - primaryRefEnd) + overlap;
        if(intronLength < mConfig.MinIntronLength)
            return RescueAttemptResult.reject(RescueRejectReason.INTRON_TOO_SHORT);
        if(intronLength > mConfig.MaxIntronLength)
            return RescueAttemptResult.reject(RescueRejectReason.INTRON_TOO_LONG);

        // Original anchors must already meet threshold — clipping only reduces them.
        final int primaryTrailingMOrig = trailingMatchedRun(primaryCigar);
        final int suppLeadingMOrig = leadingMatchedRun(suppCigar);
        if(primaryTrailingMOrig < mConfig.MinAnchorOverhang
                || suppLeadingMOrig < mConfig.MinAnchorOverhang)
            return RescueAttemptResult.reject(RescueRejectReason.SHORT_ANCHOR);

        // pick the highest-tier candidate split point in [suppLeadingS, primaryMatched]; within a
        // tier, largest min(primaryAnchor, suppAnchor) wins; ties go to trust-primary (the first
        // L seen in descending iteration). Tier-0 candidates fall through to the mate-hint and
        // trust-primary fallbacks below.
        int chosenL = -1;
        ChrIntron chosenIntron = null;
        int chosenTier = SpliceMotif.TIER_NONE;
        int chosenMinAnchor = -1;
        for(int L = primaryMatched; L >= suppLeadingS; --L)
        {
            final int primaryLoss = primaryMatched - L;
            final int suppLoss = L - suppLeadingS;
            if(primaryTrailingMOrig - primaryLoss < mConfig.MinAnchorOverhang)
                continue;
            if(suppLeadingMOrig - suppLoss < mConfig.MinAnchorOverhang)
                continue;
            if(primaryTrailingMOrig < primaryLoss || suppLeadingMOrig < suppLoss)
                continue;
            final int candidateIntronStart = primaryRefEnd - primaryLoss + 1;
            final int candidateIntronEnd = supp.Start + suppLoss - 1;
            final ChrIntron candidateIntron = new ChrIntron(
                    candidate.Chromosome, candidateIntronStart, candidateIntronEnd);
            final int tier = classifyJunctionTier(candidateIntron);
            if(tier == SpliceMotif.TIER_NONE)
                continue;
            final int candidateMinAnchor = Math.min(
                    primaryTrailingMOrig - primaryLoss,
                    suppLeadingMOrig - suppLoss);
            if(tier > chosenTier
                    || (tier == chosenTier && candidateMinAnchor > chosenMinAnchor))
            {
                chosenL = L;
                chosenIntron = candidateIntron;
                chosenTier = tier;
                chosenMinAnchor = candidateMinAnchor;
            }
        }
        // Mate-hint fallback: no annotated match, but the partner mate's rescue introduced an
        // intron whose start position falls within our overlap window. Use it as the snap target
        // instead of trust-primary fallback.
        if(chosenL == -1 && !candidate.MateHintIntrons.isEmpty())
        {
            for(ChrIntron hint : candidate.MateHintIntrons)
            {
                if(!hint.Chromosome.equals(candidate.Chromosome))
                    continue;
                final int requiredPrimaryLoss = primaryRefEnd - hint.IntronStart + 1;
                final int L = primaryMatched - requiredPrimaryLoss;
                if(L < suppLeadingS || L > primaryMatched)
                    continue;
                final int suppLoss = L - suppLeadingS;
                final int primaryLoss = requiredPrimaryLoss;
                if(primaryTrailingMOrig - primaryLoss < mConfig.MinAnchorOverhang)
                    continue;
                if(suppLeadingMOrig - suppLoss < mConfig.MinAnchorOverhang)
                    continue;
                if(primaryTrailingMOrig < primaryLoss || suppLeadingMOrig < suppLoss)
                    continue;
                final ChrIntron candidateIntron = new ChrIntron(
                        candidate.Chromosome, hint.IntronStart, supp.Start + suppLoss - 1);
                chosenL = L;
                chosenIntron = candidateIntron;
                break;
            }
        }
        if(chosenL == -1)
        {
            if(mConfig.AnnotatedOnly)
                return RescueAttemptResult.reject(RescueRejectReason.NOVEL_JUNCTION);
            // Fallback: trust primary entirely (L = primaryMatched, no clipping).
            chosenL = primaryMatched;
            chosenIntron = new ChrIntron(candidate.Chromosome, primaryRefEnd + 1, supp.Start - 1);
            // Re-check that this fallback split passes the anchor and clip constraints — supp
            // may need to clip overlap bases off its leading M.
            final int suppLoss = chosenL - suppLeadingS;
            if(suppLeadingMOrig - suppLoss < mConfig.MinAnchorOverhang || suppLeadingMOrig < suppLoss)
                return RescueAttemptResult.reject(RescueRejectReason.SHORT_ANCHOR);
        }

        final int primaryLoss = primaryMatched - chosenL;
        final int suppLoss = chosenL - suppLeadingS;

        // Build merged: primary's ops up to (excluding) trailing S, with the trailing M reduced
        // by primaryLoss; then N(intronLength); then supp's ops past the leading S, with the
        // first M reduced by suppLoss.
        final List<CigarShape.Element> merged = new ArrayList<>(primaryCigar.size() + suppCigar.size());
        for(int i = 0; i < primaryCigar.size() - 1; ++i)
        {
            if(i == primaryCigar.size() - 2 && primaryLoss > 0)
                merged.add(new CigarShape.Element(primaryCigar.get(i).Length - primaryLoss,
                        primaryCigar.get(i).Op));
            else
                merged.add(primaryCigar.get(i));
        }
        merged.add(new CigarShape.Element(intronLength, OP_SKIPPED));
        for(int i = 1; i < suppCigar.size(); ++i)
        {
            if(i == 1 && suppLoss > 0)
                merged.add(new CigarShape.Element(suppCigar.get(i).Length - suppLoss,
                        suppCigar.get(i).Op));
            else
                merged.add(suppCigar.get(i));
        }

        return RescueAttemptResult.success(primaryStart, merged, chosenIntron);
    }

    private RescueAttemptResult tryLeftExtend(
            final RescueCandidate candidate, final int primaryStart,
            final List<CigarShape.Element> primaryCigar, final RescueSupplementary supp,
            final List<CigarShape.Element> suppCigar)
    {
        final int primaryReadCovered = CigarShape.readLength(primaryCigar);
        final int suppReadCovered = CigarShape.readLength(suppCigar);

        if(primaryReadCovered != candidate.ReadLength || suppReadCovered != candidate.ReadLength)
            return RescueAttemptResult.reject(RescueRejectReason.COMPLEX_CIGAR_SHAPE);

        if(opAdjacentToSoftClip(primaryCigar, true))
            return RescueAttemptResult.reject(RescueRejectReason.COMPLEX_CIGAR_SHAPE);
        if(opAdjacentToSoftClip(suppCigar, false))
            return RescueAttemptResult.reject(RescueRejectReason.COMPLEX_CIGAR_SHAPE);

        final int primaryLeadingS = CigarShape.leadingSoftClip(primaryCigar);
        final int suppMatched = candidate.ReadLength - CigarShape.trailingSoftClip(suppCigar);
        final int overlap = suppMatched - primaryLeadingS;

        if(overlap < 0)
            return RescueAttemptResult.reject(RescueRejectReason.READ_COVERAGE_GAP);
        if(overlap > mConfig.SoftclipTolerance)
            return RescueAttemptResult.reject(RescueRejectReason.READ_COVERAGE_OVERLAP);

        final int suppRefEnd = supp.Start + CigarShape.referenceSpan(suppCigar) - 1;
        if(primaryStart <= suppRefEnd)
            return RescueAttemptResult.reject(RescueRejectReason.READ_COVERAGE_OVERLAP);

        final int intronLength = (primaryStart - 1 - suppRefEnd) + overlap;
        if(intronLength < mConfig.MinIntronLength)
            return RescueAttemptResult.reject(RescueRejectReason.INTRON_TOO_SHORT);
        if(intronLength > mConfig.MaxIntronLength)
            return RescueAttemptResult.reject(RescueRejectReason.INTRON_TOO_LONG);

        final int suppTrailingMOrig = trailingMatchedRun(suppCigar);
        final int primaryLeadingMOrig = leadingMatchedRun(primaryCigar);
        if(suppTrailingMOrig < mConfig.MinAnchorOverhang
                || primaryLeadingMOrig < mConfig.MinAnchorOverhang)
            return RescueAttemptResult.reject(RescueRejectReason.SHORT_ANCHOR);

        // mirror of tryRightExtend's snap loop. Tier-0 candidates fall through to the fallbacks below.
        int chosenL = -1;
        ChrIntron chosenIntron = null;
        int chosenTier = SpliceMotif.TIER_NONE;
        int chosenMinAnchor = -1;
        for(int L = primaryLeadingS; L <= suppMatched; ++L)
        {
            final int suppLoss = suppMatched - L;
            final int primaryLoss = L - primaryLeadingS;
            if(suppTrailingMOrig - suppLoss < mConfig.MinAnchorOverhang)
                continue;
            if(primaryLeadingMOrig - primaryLoss < mConfig.MinAnchorOverhang)
                continue;
            if(suppTrailingMOrig < suppLoss || primaryLeadingMOrig < primaryLoss)
                continue;
            final int candidateIntronStart = suppRefEnd - suppLoss + 1;
            final int candidateIntronEnd = primaryStart + primaryLoss - 1;
            final ChrIntron candidateIntron = new ChrIntron(
                    candidate.Chromosome, candidateIntronStart, candidateIntronEnd);
            final int tier = classifyJunctionTier(candidateIntron);
            if(tier == SpliceMotif.TIER_NONE)
                continue;
            final int candidateMinAnchor = Math.min(
                    suppTrailingMOrig - suppLoss,
                    primaryLeadingMOrig - primaryLoss);
            if(tier > chosenTier
                    || (tier == chosenTier && candidateMinAnchor > chosenMinAnchor))
            {
                chosenL = L;
                chosenIntron = candidateIntron;
                chosenTier = tier;
                chosenMinAnchor = candidateMinAnchor;
            }
        }
        // Mate-hint fallback for left-extend: derive L from the hinted intron's end position.
        if(chosenL == -1 && !candidate.MateHintIntrons.isEmpty())
        {
            for(ChrIntron hint : candidate.MateHintIntrons)
            {
                if(!hint.Chromosome.equals(candidate.Chromosome))
                    continue;
                final int requiredPrimaryLoss = hint.IntronEnd - (primaryStart - 1);
                final int L = primaryLeadingS + requiredPrimaryLoss;
                if(L < primaryLeadingS || L > suppMatched)
                    continue;
                final int suppLoss = suppMatched - L;
                final int primaryLoss = requiredPrimaryLoss;
                if(suppTrailingMOrig - suppLoss < mConfig.MinAnchorOverhang)
                    continue;
                if(primaryLeadingMOrig - primaryLoss < mConfig.MinAnchorOverhang)
                    continue;
                if(suppTrailingMOrig < suppLoss || primaryLeadingMOrig < primaryLoss)
                    continue;
                final ChrIntron candidateIntron = new ChrIntron(
                        candidate.Chromosome, suppRefEnd - suppLoss + 1, hint.IntronEnd);
                chosenL = L;
                chosenIntron = candidateIntron;
                break;
            }
        }
        if(chosenL == -1)
        {
            if(mConfig.AnnotatedOnly)
                return RescueAttemptResult.reject(RescueRejectReason.NOVEL_JUNCTION);
            chosenL = primaryLeadingS;
            chosenIntron = new ChrIntron(candidate.Chromosome, suppRefEnd + 1, primaryStart - 1);
            final int suppLoss = suppMatched - chosenL;
            if(suppTrailingMOrig - suppLoss < mConfig.MinAnchorOverhang || suppTrailingMOrig < suppLoss)
                return RescueAttemptResult.reject(RescueRejectReason.SHORT_ANCHOR);
        }

        final int suppLoss = suppMatched - chosenL;
        final int primaryLoss = chosenL - primaryLeadingS;

        final List<CigarShape.Element> merged = new ArrayList<>(primaryCigar.size() + suppCigar.size());
        for(int i = 0; i < suppCigar.size() - 1; ++i)
        {
            if(i == suppCigar.size() - 2 && suppLoss > 0)
                merged.add(new CigarShape.Element(suppCigar.get(i).Length - suppLoss,
                        suppCigar.get(i).Op));
            else
                merged.add(suppCigar.get(i));
        }
        merged.add(new CigarShape.Element(intronLength, OP_SKIPPED));
        for(int i = 1; i < primaryCigar.size(); ++i)
        {
            if(i == 1 && primaryLoss > 0)
                merged.add(new CigarShape.Element(primaryCigar.get(i).Length - primaryLoss,
                        primaryCigar.get(i).Op));
            else
                merged.add(primaryCigar.get(i));
        }

        return RescueAttemptResult.success(supp.Start, merged, chosenIntron);
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

    // small helpers wrapping success/failure plumbing.
    private static final class RescueAttemptResult
    {
        final RescueRejectReason Reject;
        final int MergedStart;
        final List<CigarShape.Element> MergedCigar;
        final ChrIntron IntroducedIntron;

        private RescueAttemptResult(
                final RescueRejectReason reject, final int mergedStart,
                final List<CigarShape.Element> mergedCigar, final ChrIntron intron)
        {
            Reject = reject;
            MergedStart = mergedStart;
            MergedCigar = mergedCigar;
            IntroducedIntron = intron;
        }

        static RescueAttemptResult reject(final RescueRejectReason reason)
        {
            return new RescueAttemptResult(reason, -1, null, null);
        }

        static RescueAttemptResult success(
                final int start, final List<CigarShape.Element> cigar, final ChrIntron intron)
        {
            return new RescueAttemptResult(null, start, cigar, intron);
        }
    }

    private static final class MergeAttempt
    {
        final RescueAttemptResult Inner;
        final RescueSupplementary MergedSupp;
        final boolean Ambiguous;
        final RescueRejectReason RejectReason;

        private MergeAttempt(
                final RescueAttemptResult inner, final RescueSupplementary supp,
                final boolean ambiguous, final RescueRejectReason rejectReason)
        {
            Inner = inner;
            MergedSupp = supp;
            Ambiguous = ambiguous;
            RejectReason = rejectReason;
        }

        MergeAttempt(final RescueAttemptResult inner, final RescueSupplementary supp)
        {
            this(inner, supp, false, null);
        }

        static MergeAttempt ambiguous()
        {
            return new MergeAttempt(null, null, true, RescueRejectReason.AMBIGUOUS_SUPP_CHOICE);
        }

        static MergeAttempt fail(final RescueRejectReason reason)
        {
            return new MergeAttempt(null, null, false, reason);
        }

        boolean success()
        {
            return MergedSupp != null && !Ambiguous;
        }
    }
}
