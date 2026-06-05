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

        // Supp-merge gate: a primary whose own placement is low-confidence isn't a trustworthy anchor
        // for building a spliced read from its supplementaries.
        if(candidate.PrimaryMapq < mConfig.MinPrimaryMapq)
        {
            mStatistics.countReject(RescueRejectReason.LOW_PRIMARY_MAPQ);
            return RescueResult.noMerge(RescueRejectReason.LOW_PRIMARY_MAPQ);
        }

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

        // A read can be soft-clipped at both ends (a short 5' quality clip plus a junction tail). Try
        // each clipped end as the rescue side, longer clip first (more likely the spliced tail), and
        // keep the other clip untouched in the merged cigar. First side that merges wins.
        final boolean[] sides;
        if(trailingS && leadingS)
        {
            final boolean trailingFirst =
                    CigarShape.trailingSoftClip(primaryCigar) >= CigarShape.leadingSoftClip(primaryCigar);
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
            if(sideResult.Merged)
                return sideResult;
            if(sideResult.RejectReason != RescueRejectReason.REF_VERIFY_NO_CANDIDATE_EXON)
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

    // Ref-verify one terminal softclip (rightExtend = trailing) against annotated donors/acceptors,
    // with the over-extension boundary snap. Does NOT count the no-candidate reject — tryRefVerifyOnly
    // aggregates across both ends and counts once.
    private RescueResult attemptRefVerifySide(
            final RescueCandidate candidate, final List<CigarShape.Element> primaryCigar, final boolean rightExtend)
    {
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

        // BWA often over-extends a few matched bases past the true exon boundary into the intron, so
        // the exact boundary probe misses the annotated donor/acceptor. Snap back up to MaxBoundaryShift
        // bases (smallest shift first = closest to bwa's call), trimming the over-extension into the
        // softclip, and ref-verify the enlarged softclip against the candidate exon.
        final int primaryRefEnd = candidate.PrimaryStart + CigarShape.referenceSpan(primaryCigar) - 1;
        boolean anyCandidate = false;
        RescueResult lastFailure = RescueResult.noMerge(RescueRejectReason.REF_VERIFY_MISMATCH_TOO_HIGH);

        for(int shift = 0; shift <= mConfig.MaxBoundaryShift; ++shift)
        {
            if(primaryAnchor - shift < mConfig.MinAnchorOverhang)
                break;
            final List<CigarShape.Element> shiftedCigar = shift == 0
                    ? primaryCigar
                    : shiftBoundaryIntoSoftclip(primaryCigar, shift, rightExtend);
            if(shiftedCigar == null)
                break;

            final int boundary = rightExtend ? (primaryRefEnd + 1 - shift) : (candidate.PrimaryStart - 1 + shift);
            final List<ChrIntron> candidates = rightExtend
                    ? mAnnotatedIndex.introByStart(candidate.Chromosome, boundary)
                    : mAnnotatedIndex.introByEnd(candidate.Chromosome, boundary);
            if(candidates.isEmpty())
                continue;

            anyCandidate = true;
            final RescueResult result = verifyAgainstCandidates(
                    candidate, shiftedCigar, candidates, softclipLen + shift, rightExtend);
            if(result.Merged)
                return result;
            lastFailure = result;
        }

        if(!anyCandidate)
            return RescueResult.noMerge(RescueRejectReason.REF_VERIFY_NO_CANDIDATE_EXON);
        return lastFailure;
    }

    // Moves 'shift' matched bases from the exon-edge M run into the adjacent terminal softclip,
    // realigning the boundary back toward the true exon edge. Returns null when the edge op isn't a
    // plain M run or trimming would leave nothing matched.
    private static List<CigarShape.Element> shiftBoundaryIntoSoftclip(
            final List<CigarShape.Element> cigar, final int shift, final boolean rightExtend)
    {
        final List<CigarShape.Element> shifted = new ArrayList<>(cigar);
        final int last = shifted.size() - 1;
        final int softIdx = rightExtend ? last : 0;
        final int matchIdx = rightExtend ? last - 1 : 1;
        if(matchIdx < 0 || matchIdx >= shifted.size())
            return null;

        final CigarShape.Element softElement = shifted.get(softIdx);
        final CigarShape.Element matchElement = shifted.get(matchIdx);
        if(matchElement.Op != OP_MATCH && matchElement.Op != OP_SEQ_MATCH)
            return null;
        if(matchElement.Length - shift < 1)
            return null;

        shifted.set(matchIdx, new CigarShape.Element(matchElement.Length - shift, matchElement.Op));
        shifted.set(softIdx, new CigarShape.Element(softElement.Length + shift, OP_SOFTCLIP));
        return shifted;
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

        // For each candidate exon across the junction, find the longest run of softclip bases matching
        // the exon from the JUNCTION-PROXIMAL end (within a cumulative 10% mismatch floor). bwa's softclip
        // frequently carries an outer adapter / low-quality tail that is NOT exon sequence, so only the
        // junction-proximal run is converted to M and the outer residual stays soft-clipped (e.g.
        // 19S132M -> 6S15M..N..130M, not 21M..N..130M which would force-align the adapter). The proximal
        // end is the high index of a leading softclip (the base just before the intron) and the low index
        // of a trailing softclip (just after it).
        final boolean proximalAtEnd = !rightExtend;
        ChrIntron chosen = null;
        int chosenRun = 0;
        int chosenMismatches = Integer.MAX_VALUE;
        boolean ambiguous = false;

        for(ChrIntron candidateIntron : candidates)
        {
            final int intronLength = candidateIntron.IntronEnd - candidateIntron.IntronStart + 1;
            if(intronLength < mConfig.MinIntronLength || intronLength > mConfig.MaxIntronLength)
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

            final int run = longestProximalMatchRun(softclipBases, refBases, softclipLen, proximalAtEnd);
            if(run == 0)
                continue;
            final int mismatches = mismatchesInRun(softclipBases, refBases, softclipLen, run, proximalAtEnd);

            // Longest matched run wins; tie-break on fewest mismatches. Two distinct introns tying on
            // both make the junction ambiguous (the read could splice to either) -> reject.
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
        // A partial match leaves an unexplained outer residual, so it needs a longer junction-proximal
        // anchor than a clean full match (which keeps the lenient MinAnchorOverhang floor gated upstream)
        // before we trust the spliced alignment over bwa's soft-clipped one.
        final boolean fullMatch = chosenRun == softclipLen;
        if(!fullMatch && chosenRun < mConfig.MinPartialMatchRun)
        {
            mStatistics.countReject(RescueRejectReason.REF_VERIFY_SHORT_PARTIAL_RUN);
            return RescueResult.noMerge(RescueRejectReason.REF_VERIFY_SHORT_PARTIAL_RUN);
        }

        return buildRefVerifyMerge(candidate, primaryCigar, chosen, chosenRun, softclipLen, rightExtend);
    }

    // Builds the spliced cigar from a ref-verified junction: the matched exon-proximal run becomes M, the
    // intron an N, and any outer residual of the softclip stays soft-clipped.
    private RescueResult buildRefVerifyMerge(
            final RescueCandidate candidate, final List<CigarShape.Element> primaryCigar,
            final ChrIntron chosen, final int matchedRun, final int softclipLen, final boolean rightExtend)
    {
        final int intronLength = chosen.IntronEnd - chosen.IntronStart + 1;
        final int residualSoftclip = softclipLen - matchedRun;
        final List<CigarShape.Element> merged = new ArrayList<>(primaryCigar.size() + 3);
        if(rightExtend)
        {
            for(int i = 0; i < primaryCigar.size() - 1; ++i)
                merged.add(primaryCigar.get(i));
            merged.add(new CigarShape.Element(intronLength, OP_SKIPPED));
            merged.add(new CigarShape.Element(matchedRun, OP_MATCH));
            if(residualSoftclip > 0)
                merged.add(new CigarShape.Element(residualSoftclip, OP_SOFTCLIP));
        }
        else
        {
            if(residualSoftclip > 0)
                merged.add(new CigarShape.Element(residualSoftclip, OP_SOFTCLIP));
            merged.add(new CigarShape.Element(matchedRun, OP_MATCH));
            merged.add(new CigarShape.Element(intronLength, OP_SKIPPED));
            for(int i = 1; i < primaryCigar.size(); ++i)
                merged.add(primaryCigar.get(i));
        }

        final int mergedStart = rightExtend ? candidate.PrimaryStart : (chosen.IntronStart - matchedRun);
        mStatistics.countMergedChain(1);
        return new RescueResult(
                true, CigarShape.format(merged), mergedStart,
                Collections.emptyList(), Collections.singletonList(chosen),
                1, null);
    }

    // Longest run of softclip bases matching ref from the junction-proximal end, within a cumulative
    // 10% mismatch floor (so the bases right at the splice site must be clean). The run is extended only
    // on a matching base, so it ends on a match rather than leaking a tolerated mismatch into the outer
    // residual. proximalAtEnd=true walks inward from the softclip's high index, false from its low index;
    // ref is aligned to the softclip (refBases[i] pairs with softclipBases[i]).
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
        // count viable supps per terminal side: a both-sides-clipped primary legitimately has one supp
        // for each end, so reach is ambiguous only when one side draws more than one supp.
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

        // Only merge when at most one supp is within reach of each boundary. More than one on a side
        // means the splice destination is ambiguous and we refuse to guess.
        if(rightReach > 1 || leftReach > 1)
        {
            mStatistics.countReject(RescueRejectReason.MULTIPLE_SUPPS_IN_REACH);
            return MergeOutcome.reject(RescueRejectReason.MULTIPLE_SUPPS_IN_REACH);
        }

        return best;
    }

    // Candidate ranking when (rarely) more than one supp passes — higher MAPQ wins, then smaller
    // intron. Only used to keep the best for diagnostics; >1 within reach is rejected by the caller.
    private static boolean isBetterCandidate(final MergeOutcome a, final MergeOutcome b)
    {
        if(a.MergedSupp.Mapq != b.MergedSupp.Mapq)
            return a.MergedSupp.Mapq > b.MergedSupp.Mapq;
        return intronLength(a.IntroducedIntron) < intronLength(b.IntroducedIntron);
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

        List<CigarShape.Element> suppCigar = CigarShape.parse(supp.Cigar);
        if(CigarShape.hasHardClip(suppCigar))
            return MergeOutcome.reject(RescueRejectReason.COMPLEX_CIGAR_SHAPE);

        // A tx-contig-derived supp can end up with an internal N when ContigTranslator expands a
        // cross-exon M anchor into M-N-M on the reference. If the post-N M lands inside or past the
        // primary's mapped interval, those read bases are already explained by the primary and the
        // post-N piece is a coincidental ref match. Clamp the supp to its primary-distal M anchor
        // before the rest of the merge logic looks at softclip lengths.
        final int primaryRefEnd = primaryStart + CigarShape.referenceSpan(primaryCigar) - 1;
        int suppStart = supp.Start;
        final ClampedSupp clamped = clampSuppToPrimaryBoundary(suppCigar, suppStart, primaryStart, primaryRefEnd);
        if(clamped != null)
        {
            suppCigar = clamped.Cigar;
            suppStart = clamped.Start;
            mStatistics.countSuppClamp();
        }

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
            final int suppRefEnd = suppStart + CigarShape.referenceSpan(suppCigar) - 1;
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
        return MergeOutcome.success(up.Start, merged, chosenIntron, supp, primaryIsUpstream);
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

    // Clamps a supp cigar so its post-N tail (or pre-N head, for a downstream supp) is replaced
    // with a soft clip whenever the M-after-N (or M-before-N) lands inside or past the primary's
    // mapped interval. Driven by the tx-contig-expansion artifact: ContigTranslator may split a
    // single M anchor across an exon boundary, producing M-N-M; if the inner M coincidentally
    // matches the reference inside the primary's span, the supp double-counts the same read bases
    // at two genomic positions. Returns null when no clamp is needed (no internal N, or the supp's
    // anchor sits entirely on one side of the primary without crossing back).
    private static ClampedSupp clampSuppToPrimaryBoundary(
            final List<CigarShape.Element> suppCigar, final int suppStart,
            final int primaryStart, final int primaryRefEnd)
    {
        boolean hasInternalN = false;
        for(CigarShape.Element e : suppCigar)
        {
            if(e.Op == OP_SKIPPED)
            {
                hasInternalN = true;
                break;
            }
        }
        if(!hasInternalN)
            return null;

        final int suppRefEnd = suppStart + CigarShape.referenceSpan(suppCigar) - 1;
        if(suppStart < primaryStart)
            return clampUpstreamSupp(suppCigar, suppStart, primaryStart);
        if(suppRefEnd > primaryRefEnd)
            return clampDownstreamSupp(suppCigar, suppStart, primaryRefEnd);
        return null;
    }

    // Walk left-to-right; remember the first N. When a subsequent M's ref start is at or past
    // primaryStart, cut the cigar at that first N (inclusive) and replace the dropped tail with a
    // trailing softclip sized to the dropped read bases.
    private static ClampedSupp clampUpstreamSupp(
            final List<CigarShape.Element> suppCigar, final int suppStart, final int primaryStart)
    {
        int refCursor = suppStart;
        int readCursor = 0;
        int firstNIdx = -1;
        int readAtFirstN = 0;

        for(int i = 0; i < suppCigar.size(); i++)
        {
            final CigarShape.Element e = suppCigar.get(i);
            final char op = e.Op;
            final int len = e.Length;

            if(op == OP_SKIPPED && firstNIdx == -1)
            {
                firstNIdx = i;
                readAtFirstN = readCursor;
            }

            if(firstNIdx != -1 && refCursor >= primaryStart
                    && (op == OP_MATCH || op == OP_SEQ_MATCH || op == OP_SEQ_MISMATCH))
            {
                return cutTailAt(suppCigar, suppStart, firstNIdx, readAtFirstN);
            }

            if(CigarShape.consumesReference(op))
                refCursor += len;
            if(CigarShape.consumesRead(op))
                readCursor += len;
        }
        return null;
    }

    // Mirror of clampUpstreamSupp: find the last N whose preceding M ended at or before
    // primaryRefEnd, drop everything up to and including that N, and add a leading softclip sized
    // to the dropped read bases. New start = the position of the first kept reference op.
    private static ClampedSupp clampDownstreamSupp(
            final List<CigarShape.Element> suppCigar, final int suppStart, final int primaryRefEnd)
    {
        int refCursor = suppStart;
        int readCursor = 0;
        int lastBadNIdx = -1;
        int readAfterLastBadN = 0;
        int refAfterLastBadN = -1;
        boolean prevMInsidePrimary = false;

        for(int i = 0; i < suppCigar.size(); i++)
        {
            final CigarShape.Element e = suppCigar.get(i);
            final char op = e.Op;
            final int len = e.Length;

            if(op == OP_MATCH || op == OP_SEQ_MATCH || op == OP_SEQ_MISMATCH)
            {
                final int mRefEnd = refCursor + len - 1;
                prevMInsidePrimary = mRefEnd <= primaryRefEnd;
            }

            if(op == OP_SKIPPED && prevMInsidePrimary)
            {
                lastBadNIdx = i;
                // ref/read cursors just after this N op
                refAfterLastBadN = refCursor + len;
                readAfterLastBadN = readCursor;
            }

            if(CigarShape.consumesReference(op))
                refCursor += len;
            if(CigarShape.consumesRead(op))
                readCursor += len;
        }

        if(lastBadNIdx == -1 || refAfterLastBadN < 0)
            return null;
        return cutHeadAt(suppCigar, refAfterLastBadN, lastBadNIdx, readAfterLastBadN);
    }

    private static ClampedSupp cutTailAt(
            final List<CigarShape.Element> suppCigar, final int suppStart,
            final int cutIdx, final int readBeforeCut)
    {
        final List<CigarShape.Element> trimmed = new ArrayList<>(cutIdx + 1);
        for(int i = 0; i < cutIdx; i++)
            trimmed.add(suppCigar.get(i));
        final int totalRead = CigarShape.readLength(suppCigar);
        final int trailingS = totalRead - readBeforeCut;
        if(trailingS > 0)
            trimmed.add(new CigarShape.Element(trailingS, OP_SOFTCLIP));
        return new ClampedSupp(suppStart, trimmed);
    }

    private static ClampedSupp cutHeadAt(
            final List<CigarShape.Element> suppCigar, final int newStart,
            final int cutIdx, final int readBeforeCut)
    {
        final List<CigarShape.Element> trimmed = new ArrayList<>(suppCigar.size() - cutIdx);
        if(readBeforeCut > 0)
            trimmed.add(new CigarShape.Element(readBeforeCut, OP_SOFTCLIP));
        for(int i = cutIdx + 1; i < suppCigar.size(); i++)
            trimmed.add(suppCigar.get(i));
        return new ClampedSupp(newStart, trimmed);
    }

    private static final class ClampedSupp
    {
        final int Start;
        final List<CigarShape.Element> Cigar;

        ClampedSupp(final int start, final List<CigarShape.Element> cigar)
        {
            Start = start;
            Cigar = cigar;
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

    // Result of attempting to merge one supp into the primary. Two shapes: success (carries the
    // merged geometry plus the supp identity) or reject (carries a reason).
    private static final class MergeOutcome
    {
        final RescueRejectReason Reject;
        final int MergedStart;
        final List<CigarShape.Element> MergedCigar;
        final ChrIntron IntroducedIntron;
        final RescueSupplementary MergedSupp;
        final boolean RightExtend;       // which terminal softclip of the primary this merge resolved

        private MergeOutcome(
                final RescueRejectReason reject,
                final int mergedStart, final List<CigarShape.Element> mergedCigar,
                final ChrIntron intron, final RescueSupplementary supp, final boolean rightExtend)
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
                final int start, final List<CigarShape.Element> cigar,
                final ChrIntron intron, final RescueSupplementary supp, final boolean rightExtend)
        {
            return new MergeOutcome(null, start, cigar, intron, supp, rightExtend);
        }

        boolean isSuccess()
        {
            return Reject == null;
        }
    }
}
