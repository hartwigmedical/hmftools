package com.hartwig.hmftools.tars.liftback.supplementary;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.tars.common.BwaScoring;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

// Merges a primary's terminal softclip with an annotated-intron-spanning supplementary into a
// single spliced primary (N op). Chains up to SupplementaryConfig.MaxSuppMerges merges per read; README Step 3.
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
    // supps. On failure, rejectReason carries the gate hit.
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
        MULTIPLE_SUPPS_IN_REACH       // more than one supp within merge reach; refuse to guess which splice
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
            return Result.noMerge(RejectReason.NO_MATCHING_SUPP);
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
        ClampedSupp clamped = clampSuppToPrimaryBoundary(suppCigar, suppStart, primaryStart, primarySide.RefEnd);
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

    // Direction-agnostic merge: validates the anchor pair, scores candidate junction positions by tier, falls
    // back to the mate hint then the midpoint. supp is threaded as result payload only.
    private MergeOutcome mergeJunction(
            final Candidate candidate, final Side up, final Side down,
            final boolean primaryIsUpstream, final Supplementary supp)
    {
        if(CigarUtils.cigarBaseLength(up.Cigar) != candidate.readLength()
                || CigarUtils.cigarBaseLength(down.Cigar) != candidate.readLength())
            return MergeOutcome.reject(RejectReason.COMPLEX_CIGAR_SHAPE);

        if(opAdjacentToSoftClip(up.Cigar, false) || opAdjacentToSoftClip(down.Cigar, true))
        {
            return MergeOutcome.reject(RejectReason.COMPLEX_CIGAR_SHAPE);
        }

        int upMatchedRead = candidate.readLength() - up.TrailingS;
        int overlap = upMatchedRead - down.LeadingS;

        if(overlap < 0)
        {
            return MergeOutcome.reject(RejectReason.READ_COVERAGE_GAP);
        }
        if(overlap > mConfig.MaxSuppReadOverlap)
        {
            return MergeOutcome.reject(RejectReason.READ_COVERAGE_OVERLAP);
        }
        if(down.Start <= up.RefEnd)
        {
            return MergeOutcome.reject(RejectReason.READ_COVERAGE_OVERLAP);
        }

        // Intron length is invariant under the junction position - check once up front.
        int intronLength = (down.Start - 1 - up.RefEnd) + overlap;
        if(intronLength < mConfig.MinIntronLength)
        {
            return MergeOutcome.reject(RejectReason.INTRON_TOO_SHORT);
        }
        if(intronLength > mConfig.MaxIntronLength)
        {
            return MergeOutcome.reject(RejectReason.INTRON_TOO_LONG);
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
        List<CigarElement> merged = buildMergedCigar(up.Cigar, down.Cigar, upLoss, downLoss, intronLength);
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

    // Length of the M/=/X run adjacent to the terminal softclip (e.g. "5M 50S" -> 5). fromEnd scans the trailing side.
    private static int matchedRun(final List<CigarElement> elements, final boolean fromEnd)
    {
        int start = fromEnd ? elements.size() - 1 : 0;
        int step = fromEnd ? -1 : 1;
        for(int i = start; i >= 0 && i < elements.size(); i += step)
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

    // ContigTranslator can expand a cross-exon M into M-N-M, leaving a supp that spans the junction itself and so
    // carries no terminal softclip to pair the primary with. Split it at an internal N and keep only the block away
    // from the primary, softclipping the rest, so it can serve as one anchor of the merge. Null when there is nothing
    // to cut: a supp starting before the primary keeps its head, one overrunning the primary's end keeps its tail.
    private static ClampedSupp clampSuppToPrimaryBoundary(
            final List<CigarElement> suppCigar, final int suppStart,
            final int primaryStart, final int primaryRefEnd)
    {
        boolean keepHead = suppStart < primaryStart;
        if(!keepHead && suppStart + CigarUtils.cigarAlignedLength(suppCigar) - 1 <= primaryRefEnd)
        {
            return null;
        }

        int refCursor = suppStart;
        int readCursor = 0;
        int splitIndex = -1;
        int readAtSplit = 0;
        int refAfterSplit = -1;
        boolean lastBlockInsidePrimary = false;

        for(int i = 0; i < suppCigar.size(); ++i)
        {
            CigarOperator op = suppCigar.get(i).getOperator();
            int length = suppCigar.get(i).getLength();

            if(keepHead)
            {
                // split at the first N, but only once a later block reaches into the primary's span
                if(op == CigarOperator.N && splitIndex == -1)
                {
                    splitIndex = i;
                    readAtSplit = readCursor;
                }
                else if(splitIndex != -1 && refCursor >= primaryStart && op.isAlignment())
                {
                    return cutAt(suppCigar, suppStart, splitIndex, readAtSplit, keepHead);
                }
            }
            else
            {
                // split at the last N whose preceding block still ended inside the primary's span
                if(op.isAlignment())
                {
                    lastBlockInsidePrimary = refCursor + length - 1 <= primaryRefEnd;
                }
                else if(op == CigarOperator.N && lastBlockInsidePrimary)
                {
                    splitIndex = i;
                    readAtSplit = readCursor;
                    refAfterSplit = refCursor + length;
                }
            }

            if(op.consumesReferenceBases())
            {
                refCursor += length;
            }
            if(op.consumesReadBases())
            {
                readCursor += length;
            }
        }

        // the keep-head cut returns from inside the loop, so reaching here with one pending means no block ever
        // reached the primary
        if(keepHead || splitIndex == -1)
        {
            return null;
        }

        return cutAt(suppCigar, refAfterSplit, splitIndex, readAtSplit, keepHead);
    }

    // Drops the split N and everything on the far side of it, replaced by a softclip over the read bases that side
    // covered. The kept side keeps its own coordinates, so only a kept tail needs a new start.
    private static ClampedSupp cutAt(
            final List<CigarElement> suppCigar, final int start,
            final int splitIndex, final int readAtSplit, final boolean keepHead)
    {
        List<CigarElement> trimmed = new ArrayList<>(suppCigar.size());

        if(keepHead)
        {
            trimmed.addAll(suppCigar.subList(0, splitIndex));
            int clipLength = CigarUtils.cigarBaseLength(suppCigar) - readAtSplit;
            if(clipLength > 0)
            {
                trimmed.add(new CigarElement(clipLength, CigarOperator.S));
            }
        }
        else
        {
            if(readAtSplit > 0)
            {
                trimmed.add(new CigarElement(readAtSplit, CigarOperator.S));
            }
            trimmed.addAll(suppCigar.subList(splitIndex + 1, suppCigar.size()));
        }

        return new ClampedSupp(start, trimmed);
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
            return new Side(
                    start, cigar,
                    CigarUtils.leftSoftClipLength(cigar), CigarUtils.rightSoftClipLength(cigar),
                    matchedRun(cigar, false), matchedRun(cigar, true),
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
