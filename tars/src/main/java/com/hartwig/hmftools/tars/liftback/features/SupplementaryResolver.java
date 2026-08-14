package com.hartwig.hmftools.tars.liftback.features;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static com.hartwig.hmftools.tars.common.TarsCigarUtils.indelAdjacentToTerminalSoftClip;
import static com.hartwig.hmftools.tars.common.TarsCigarUtils.terminalMatchedRun;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.tars.liftback.EnsemblAnnotationIndex;
import com.hartwig.hmftools.tars.common.BwaScoring;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

// Merges a primary candidate's terminal softclip with a supplementary into a spliced placement candidate.
//
// A supplementary merges into the primary only when all of these hold; RejectReason names the first that fails.
//   1. the primary has a terminal softclip to extend across      NO_TERMINAL_SOFTCLIP
//   2. exactly one supplementary reaches that side               MULTIPLE_SUPPS_IN_REACH
//   3. the supp is on the same chromosome and the same strand    DIFFERENT_CHROMOSOME, OPPOSITE_STRAND
//   4. the cigars are complementary and simple - no hard clip,
//      no indel at the boundary, full read length, one clear side  NO_MATCHING_SUPP, COMPLEX_CIGAR_SHAPE
//   5. primary M and supp M abut on the read, overlapping by at
//      most MaxSuppReadOverlap and not at all on the reference     READ_COVERAGE_GAP, READ_COVERAGE_OVERLAP
//   6. the implied intron is within [MinIntronLength, MaxIntronLength]  INTRON_TOO_SHORT, INTRON_TOO_LONG
//   7. both anchors are at least MinAnchorOverhang                 SHORT_ANCHOR
//   8. a junction position scores above Tier.NONE, or is annotated
//      when AnnotatedOnly is set                                   NOVEL_JUNCTION
public class SupplementaryResolver
{
    // readBases seeds the tie-break between equally scoring junction positions; mateHintIntrons biases that choice
    // toward the partner mate's resolved junctions.
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

    // index is the supplementary's position in the caller's list, so the resolver can report which supps to drop.
    public record Supplementary(
            int index, String chromosome, boolean forwardStrand, int start, String cigar, int mapQuality)
    {
    }

    // mergedCigar/mergedStart describe the new primary and droppedSupplementaryIndices the absorbed supps on success;
    // rejectReason carries the gate that was hit on failure.
    public record Result(
            boolean merged, String mergedCigar, int mergedStart, List<Integer> droppedSupplementaryIndices,
            List<ChrBaseRegion> introducedIntrons, int chainDepth, int spliceStrand, RejectReason rejectReason)
    {
        private static final Result NO_MERGE_NO_OP = new Result(
                false, null, -1, Collections.emptyList(), Collections.emptyList(), 0, 0, RejectReason.NO_TERMINAL_SOFTCLIP);

        public static Result noMerge(final RejectReason reason)
        {
            if(reason == RejectReason.NO_TERMINAL_SOFTCLIP)
            {
                return NO_MERGE_NO_OP;
            }
            return new Result(false, null, -1, Collections.emptyList(), Collections.emptyList(), 0, 0, reason);
        }
    }

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
    // are ranked with compareTo.
    public enum Tier
    {
        NONE,            // neither donor nor acceptor matches a known splice motif
        SEMI_CANONICAL,  // GC-AG / AT-AC (and reverse-complement equivalents)
        CANONICAL,       // GT-AG (~99% of splice sites)
        ANNOTATED        // matches an annotated junction from the sidecar
    }

    // no junction position scored above Tier.NONE; a real position is a read offset, so never negative
    private static final int NO_JUNCTION_POSITION = -1;

    private final EnsemblAnnotationIndex mEnsemblAnnotationIndex;
    private final RefGenomeInterface mRefGenome;
    private final SupplementaryConfig mConfig;

    public SupplementaryResolver(final Set<ChrBaseRegion> annotatedJunctions, final SupplementaryConfig config)
    {
        this(
                EnsemblAnnotationIndex.fromJunctions(annotatedJunctions != null ? annotatedJunctions : new HashSet<>()),
                null, config);
    }

    public SupplementaryResolver(
            final EnsemblAnnotationIndex annotationIndex, final RefGenomeInterface refGenome,
            final SupplementaryConfig config)
    {
        mEnsemblAnnotationIndex = annotationIndex != null ? annotationIndex : EnsemblAnnotationIndex.fromJunctions(new HashSet<>());
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
        int spliceStrand = 0;
        boolean conflictingStrands = false;

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
            if(merge.SpliceStrand != 0)
            {
                if(spliceStrand == 0 && !conflictingStrands)
                {
                    spliceStrand = merge.SpliceStrand;
                }
                else if(spliceStrand != merge.SpliceStrand)
                {
                    spliceStrand = 0;
                    conflictingStrands = true;
                }
            }
            remaining.remove(merge.MergedSupp);
            ++chainDepth;
        }

        if(chainDepth == 0)
        {
            // unreachable: a first-iteration failure returns above
            return Result.noMerge(lastReject != null ? lastReject : RejectReason.NO_MATCHING_SUPP);
        }

        return new Result(
                true, CigarUtils.cigarElementsToStr(primaryCigar), primaryStart,
                dropped, introns, chainDepth, spliceStrand, null);
    }

    private Tier classifyJunctionTier(final ChrBaseRegion candidateIntron)
    {
        if(mEnsemblAnnotationIndex.containsJunction(candidateIntron))
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

    private int junctionStrand(final ChrBaseRegion intron)
    {
        int annotatedStrand = mEnsemblAnnotationIndex.junctionStrand(intron);
        if(annotatedStrand != 0 || mRefGenome == null)
        {
            return annotatedStrand;
        }

        byte[] donor = refBases(intron.Chromosome, intron.start(), intron.start() + 1);
        byte[] acceptor = refBases(intron.Chromosome, intron.end() - 1, intron.end());
        return motifStrand(donor, acceptor);
    }

    public int spliceStrand(final String chromosome, final int start, final String cigar)
    {
        int referencePosition = start;
        int strand = 0;
        for(CigarElement element : CigarUtils.cigarElementsFromStr(cigar))
        {
            if(element.getOperator() == CigarOperator.N)
            {
                int intronStrand = junctionStrand(new ChrBaseRegion(
                        chromosome, referencePosition, referencePosition + element.getLength() - 1));
                if(intronStrand == 0 || strand != 0 && strand != intronStrand)
                {
                    return 0;
                }
                strand = intronStrand;
            }
            if(element.getOperator().consumesReferenceBases())
            {
                referencePosition += element.getLength();
            }
        }
        return strand;
    }

    static int motifStrand(final byte[] donorBases, final byte[] acceptorBases)
    {
        if(donorBases == null || donorBases.length != 2 || acceptorBases == null || acceptorBases.length != 2)
        {
            return 0;
        }

        String donor = upperCase(donorBases);
        String acceptor = upperCase(acceptorBases);
        if(donor.equals("GT") && acceptor.equals("AG")
                || donor.equals("GC") && acceptor.equals("AG")
                || donor.equals("AT") && acceptor.equals("AC"))
        {
            return 1;
        }
        if(donor.equals("CT") && acceptor.equals("AC")
                || donor.equals("CT") && acceptor.equals("GC")
                || donor.equals("GT") && acceptor.equals("AT"))
        {
            return -1;
        }
        return 0;
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

            // At most one supp per terminal softclip can be resolved. A second reaching the same side leaves no way to
            // tell which splice is real, so reject the read rather than guess.
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

    // Higher MAPQ wins, then the smaller intron. Only ever compares a right-extend against a left-extend: a second supp
    // on either side is rejected above, and the chain loop picks up the loser next pass.
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

        // ContigTranslator can expand a cross-exon M into M-N-M. When the post-N M overlaps the primary's span, clamp
        // the supp to its primary-distal anchor first.
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

        // Clipped at both ends, so the supp could extend either way. The side of the primary it sits on settles it; if
        // it sits cleanly on neither, refuse to guess.
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

    // Conditions 4-6 of the merge policy for a resolved up/down anchor pair, cheapest first. Returns the first failure,
    // or null when the pair may merge.
    private RejectReason anchorPairReject(
            final Candidate candidate, final Side up, final Side down, final int overlap, final int intronLength)
    {
        if(CigarUtils.cigarBaseLength(up.Cigar) != candidate.readLength()
                || CigarUtils.cigarBaseLength(down.Cigar) != candidate.readLength())
        {
            return RejectReason.COMPLEX_CIGAR_SHAPE;
        }

        if(indelAdjacentToTerminalSoftClip(up.Cigar, false) || indelAdjacentToTerminalSoftClip(down.Cigar, true))
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

    // Direction-agnostic merge: validates the anchor pair, scores junction positions by tier, falls back to the mate
    // hint then the midpoint. supp is carried through as result payload only.
    private MergeOutcome mergeJunction(
            final Candidate candidate, final Side up, final Side down,
            final boolean primaryIsUpstream, final Supplementary supp)
    {
        int upMatchedRead = candidate.readLength() - up.TrailingS;
        int overlap = upMatchedRead - down.LeadingS;

        // Intron length is invariant under the junction position, so compute it once.
        int intronLength = (down.Start - 1 - up.RefEnd) + overlap;

        RejectReason gateReject = anchorPairReject(candidate, up, down, overlap, intronLength);
        if(gateReject != null)
        {
            return MergeOutcome.reject(gateReject);
        }

        // Junction position priority: annotated boundary (sidecar) > canonical > semi-canonical motif, then the mate's
        // intron, then the midpoint of the ambiguous range.
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
            // paths already check the anchors, so only this fallback has to.
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
        ChrBaseRegion intron = intronAt(candidate, up, down, upMatchedRead, junctionReadPosition);
        return MergeOutcome.success(up.Start, merged, intron, supp, primaryIsUpstream, junctionStrand(intron));
    }

    // Intron implied by putting the junction at this read position: the upstream alignment gives up the bases past it,
    // the downstream alignment the bases before it.
    private static ChrBaseRegion intronAt(
            final Candidate candidate, final Side up, final Side down,
            final int upMatchedRead, final int readPosition)
    {
        int upLoss = upMatchedRead - readPosition;
        int downLoss = readPosition - down.LeadingS;
        return new ChrBaseRegion(
                candidate.chromosome(), up.RefEnd - upLoss + 1, down.Start + downLoss - 1);
    }

    // Scores every valid junction position by tier and returns one at the highest tier present. Ties there are broken
    // pseudo-randomly but deterministically (read-seeded) so ambiguous junctions distribute yet stay reproducible.
    // Returns NO_JUNCTION_POSITION when no position reaches a splice motif or annotated boundary.
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

    // Deterministic per-read seed: hashing the read bases (when available) with the upstream boundary keeps a tie
    // reproducible run to run while letting two reads over the same junction pick differently.
    private static int tieBreakSeed(final Candidate candidate, final Side up)
    {
        int base = candidate.readBases() != null ? Arrays.hashCode(candidate.readBases()) : 0;
        return 31 * base + up.RefEnd;
    }

    // Use the mate's resolved intron as a hint: primary-upstream pins the intron start, primary-downstream pins its
    // end. Either way the pinned end fixes the read position, and intronAt reproduces the hinted intron.
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
        for(int i = 0; i < upCigar.size() - 1; ++i)
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
        for(int i = 1; i < downCigar.size(); ++i)
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
                if(op == CigarOperator.N && splitIndex == -1)
                {
                    splitIndex = i;
                    readAtSplit = readCursor;
                }
                else if(splitIndex != -1 && refCursor >= primaryStart && op.isAlignment())
                {
                    return cutAt(suppCigar, suppStart, splitIndex, readAtSplit, true);
                }
            }
            else
            {
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

        if(keepHead || splitIndex == -1)
        {
            return null;
        }

        return cutAt(suppCigar, refAfterSplit, splitIndex, readAtSplit, false);
    }

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
                    terminalMatchedRun(cigar, false), terminalMatchedRun(cigar, true),
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
        int SpliceStrand;

        private MergeOutcome(
                final RejectReason reject,
                final int mergedStart, final List<CigarElement> mergedCigar,
                final ChrBaseRegion intron, final Supplementary supp, final boolean rightExtend, final int spliceStrand)
        {
            Reject = reject;
            MergedStart = mergedStart;
            MergedCigar = mergedCigar;
            IntroducedIntron = intron;
            MergedSupp = supp;
            RightExtend = rightExtend;
            SpliceStrand = spliceStrand;
        }

        static MergeOutcome reject(final RejectReason reason)
        {
            return new MergeOutcome(reason, -1, null, null, null, false, 0);
        }

        static MergeOutcome success(
                final int start, final List<CigarElement> cigar,
                final ChrBaseRegion intron, final Supplementary supp, final boolean rightExtend, final int spliceStrand)
        {
            return new MergeOutcome(null, start, cigar, intron, supp, rightExtend, spliceStrand);
        }

        boolean isSuccess()
        {
            return Reject == null;
        }
    }
}
