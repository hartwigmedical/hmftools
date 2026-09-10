package com.hartwig.hmftools.tars.liftback.features;

import static com.hartwig.hmftools.tars.common.TarsCigarUtils.indelAdjacentToTerminalSoftClip;
import static com.hartwig.hmftools.tars.common.TarsCigarUtils.terminalMatchedRun;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.liftback.EnsemblAnnotationIndex;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class SupplementaryMerger
{
    // NOTE: generated equals/hashCode compare readBases by array identity, so do not use as a Set/Map key.
    public record Placement(
            String chromosome, boolean forwardStrand, int readLength, int primaryStart, String primaryCigar,
            List<Supplementary> supplementaries, byte[] readBases, List<ChrBaseRegion> mateHintIntrons)
    {
        public Placement
        {
            if(mateHintIntrons == null)
            {
                mateHintIntrons = Collections.emptyList();
            }
        }

        public Placement(
                final String chromosome, final boolean forwardStrand, final int readLength,
                final int primaryStart, final String primaryCigar, final List<Supplementary> supplementaries)
        {
            this(chromosome, forwardStrand, readLength, primaryStart, primaryCigar, supplementaries, null, Collections.emptyList());
        }
    }

    public record Supplementary(
            int index, String chromosome, boolean forwardStrand, int start, String cigar, int mapQuality)
    {
    }

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
        NO_TERMINAL_SOFTCLIP,
        NO_MATCHING_SUPP,
        DIFFERENT_CHROMOSOME,
        OPPOSITE_STRAND,
        READ_COVERAGE_OVERLAP,
        READ_COVERAGE_GAP,
        INTRON_TOO_SHORT,
        INTRON_TOO_LONG,
        SHORT_ANCHOR,
        NOVEL_JUNCTION,
        COMPLEX_CIGAR_SHAPE,
        MULTIPLE_SUPPS_IN_REACH
    }

    private static final int NO_JUNCTION_POSITION = -1;

    private final SpliceJunctions mSpliceJunctions;
    private final EnsemblAnnotationIndex mEnsemblAnnotationIndex;
    private final SupplementaryConfig mConfig;

    public SupplementaryMerger(final Set<ChrBaseRegion> annotatedJunctions, final SupplementaryConfig config)
    {
        this(
                EnsemblAnnotationIndex.fromJunctions(annotatedJunctions != null ? annotatedJunctions : new HashSet<>()),
                null, config);
    }

    public SupplementaryMerger(
            final EnsemblAnnotationIndex annotationIndex, final RefGenomeInterface refGenome,
            final SupplementaryConfig config)
    {
        mEnsemblAnnotationIndex = annotationIndex != null
                ? annotationIndex : EnsemblAnnotationIndex.fromJunctions(new HashSet<>());
        mSpliceJunctions = new SpliceJunctions(mEnsemblAnnotationIndex, refGenome);
        mConfig = config;
    }

    public int spliceStrand(final String chromosome, final int start, final String cigar)
    {
        return mSpliceJunctions.spliceStrand(chromosome, start, cigar);
    }

    public Result merge(final Placement placement)
    {
        if(placement.supplementaries().isEmpty())
        {
            return Result.noMerge(RejectReason.NO_MATCHING_SUPP);
        }

        List<CigarElement> primaryCigar = CigarUtils.cigarElementsFromStr(placement.primaryCigar());
        if(CigarUtils.hasHardClip(primaryCigar))
        {
            return Result.noMerge(RejectReason.COMPLEX_CIGAR_SHAPE);
        }

        int primaryStart = placement.primaryStart();
        List<Supplementary> remaining = new ArrayList<>(placement.supplementaries());
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
                    placement, primaryStart, primaryCigar, remaining);

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
            // Do not reuse another XA placement from the absorbed supplementary record.
            remaining.removeIf(supp -> supp.index() == merge.MergedSupp.index());
            ++chainDepth;
        }

        if(chainDepth == 0)
        {
            return Result.noMerge(lastReject != null ? lastReject : RejectReason.NO_MATCHING_SUPP);
        }

        return new Result(
                true, CigarUtils.cigarElementsToStr(primaryCigar), primaryStart,
                dropped, introns, chainDepth, spliceStrand, null);
    }

    private MergeOutcome pickBestSupplementary(
            final Placement placement, final int primaryStart,
            final List<CigarElement> primaryCigar, final List<Supplementary> supps)
    {
        List<Supplementary> selectedSupplementaries =
                selectSupplementaryPlacements(placement, primaryStart, primaryCigar, supps);

        MergeOutcome chosen = null;
        RejectReason lastReject = null;
        boolean rightMerged = false;
        boolean leftMerged = false;
        for(Supplementary supp : selectedSupplementaries)
        {
            MergeOutcome outcome = tryMerge(placement, primaryStart, primaryCigar, supp);
            if(!outcome.isSuccess())
            {
                lastReject = outcome.Reject;
                continue;
            }

            // Reject multiple valid alignments reaching the same soft clip.
            if(outcome.RightExtend ? rightMerged : leftMerged)
            {
                return MergeOutcome.reject(RejectReason.MULTIPLE_SUPPS_IN_REACH);
            }

            rightMerged |= outcome.RightExtend;
            leftMerged |= !outcome.RightExtend;

            if(chosen == null || isBetterMerge(outcome, chosen))
            {
                chosen = outcome;
            }
        }

        return chosen != null
                ? chosen
                : MergeOutcome.reject(lastReject != null ? lastReject : RejectReason.NO_MATCHING_SUPP);
    }

    // Select one placement per supplementary record before attempting a merge.
    public static List<Supplementary> selectSupplementaryPlacements(final Placement placement)
    {
        return selectSupplementaryPlacements(
                placement, placement.primaryStart(), CigarUtils.cigarElementsFromStr(placement.primaryCigar()),
                placement.supplementaries());
    }

    static List<Supplementary> selectSupplementaryPlacements(
            final Placement placement, final int primaryStart, final List<CigarElement> primaryCigar,
            final List<Supplementary> supplementaries)
    {
        Map<Integer, Supplementary> selected = new LinkedHashMap<>();

        for(Supplementary candidate : supplementaries)
        {
            Supplementary current = selected.get(candidate.index());
            if(current == null)
            {
                selected.put(candidate.index(), candidate);
                continue;
            }

            if(current.mapQuality() != 0)
            {
                continue;
            }

            if(comparePlacementPriority(placement, primaryStart, primaryCigar, candidate, current) < 0)
            {
                selected.put(candidate.index(), candidate);
            }
        }

        return new ArrayList<>(selected.values());
    }

    private static int comparePlacementPriority(
            final Placement placement, final int primaryStart, final List<CigarElement> primaryCigar,
            final Supplementary candidate, final Supplementary current)
    {
        PlacementPairPriority candidatePriority = placementPriority(
                placement, primaryStart, primaryCigar, candidate);
        PlacementPairPriority currentPriority = placementPriority(
                placement, primaryStart, primaryCigar, current);

        int priorityComparison = candidatePriority.compareTo(currentPriority);
        if(priorityComparison != 0)
        {
            return priorityComparison;
        }

        int candidateRandomOrder = 31 * Arrays.hashCode(placement.readBases()) + candidate.hashCode();
        int currentRandomOrder = 31 * Arrays.hashCode(placement.readBases()) + current.hashCode();
        return Integer.compareUnsigned(candidateRandomOrder, currentRandomOrder);
    }

    private static PlacementPairPriority placementPriority(
            final Placement placement, final int primaryStart, final List<CigarElement> primaryCigar,
            final Supplementary supplementary)
    {
        if(placement.chromosome().equals(supplementary.chromosome()))
        {
            Side primarySide = Side.of(primaryStart, primaryCigar);
            Side supplementarySide = Side.of(
                    supplementary.start(), CigarUtils.cigarElementsFromStr(supplementary.cigar()));
            if(primarySide.LeadingS != supplementarySide.LeadingS)
            {
                boolean primaryLinksEnd = primarySide.LeadingS < supplementarySide.LeadingS;
                boolean supplementaryLinksEnd = !primaryLinksEnd;

                int primaryBreakend = breakendPosition(primarySide, placement.forwardStrand(), primaryLinksEnd);
                int supplementaryBreakend =
                        breakendPosition(supplementarySide, supplementary.forwardStrand(), supplementaryLinksEnd);
                Orientation primaryOrientation = breakendOrientation(placement.forwardStrand(), primaryLinksEnd);
                Orientation supplementaryOrientation =
                        breakendOrientation(supplementary.forwardStrand(), supplementaryLinksEnd);

                return PlacementPairPriority.between(
                        placement.chromosome(), primaryBreakend, primaryOrientation,
                        supplementary.chromosome(), supplementaryBreakend, supplementaryOrientation);
            }
        }
        return PlacementPairPriority.fallback();
    }

    private static int breakendPosition(final Side side, final boolean forwardStrand, final boolean linksEnd)
    {
        if(linksEnd)
        {
            return forwardStrand ? side.RefEnd : side.Start;
        }
        return forwardStrand ? side.Start : side.RefEnd;
    }

    private static Orientation breakendOrientation(final boolean forwardStrand, final boolean linksEnd)
    {
        return linksEnd == forwardStrand ? Orientation.FORWARD : Orientation.REVERSE;
    }

    // Higher MAPQ wins, then the smaller intron.
    private static boolean isBetterMerge(final MergeOutcome outcome, final MergeOutcome chosen)
    {
        if(outcome.MergedSupp.mapQuality() != chosen.MergedSupp.mapQuality())
        {
            return outcome.MergedSupp.mapQuality() > chosen.MergedSupp.mapQuality();
        }
        return outcome.IntroducedIntron.baseLength() < chosen.IntroducedIntron.baseLength();
    }

    private MergeOutcome tryMerge(
            final Placement placement, final int primaryStart,
            final List<CigarElement> primaryCigar, final Supplementary supp)
    {
        if(!placement.chromosome().equals(supp.chromosome()))
        {
            return MergeOutcome.reject(RejectReason.DIFFERENT_CHROMOSOME);
        }

        if(placement.forwardStrand() != supp.forwardStrand())
        {
            return MergeOutcome.reject(RejectReason.OPPOSITE_STRAND);
        }

        List<CigarElement> suppCigar = CigarUtils.cigarElementsFromStr(supp.cigar());
        if(CigarUtils.hasHardClip(suppCigar))
        {
            return MergeOutcome.reject(RejectReason.COMPLEX_CIGAR_SHAPE);
        }

        // Keep the primary-distal block of a translated M-N-M alignment.
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

        // Genomic position resolves a supplementary clipped at both ends.
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

        return mergeJunction(placement, up, down, primaryIsUpstream, supp);
    }

    private RejectReason anchorPairReject(
            final Placement placement, final Side up, final Side down, final int overlap, final int intronLength)
    {
        if(CigarUtils.cigarBaseLength(up.Cigar) != placement.readLength()
                || CigarUtils.cigarBaseLength(down.Cigar) != placement.readLength())
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

    private MergeOutcome mergeJunction(
            final Placement placement, final Side up, final Side down,
            final boolean primaryIsUpstream, final Supplementary supp)
    {
        int upMatchedRead = placement.readLength() - up.TrailingS;
        int overlap = upMatchedRead - down.LeadingS;

        int intronLength = (down.Start - 1 - up.RefEnd) + overlap;

        RejectReason gateReject = anchorPairReject(placement, up, down, overlap, intronLength);
        if(gateReject != null)
        {
            return MergeOutcome.reject(gateReject);
        }

        // Annotation, motif, mate hint, then midpoint.
        int junctionReadPosition = scanJunctionPositions(placement, up, down, upMatchedRead, primaryIsUpstream);

        if(junctionReadPosition == NO_JUNCTION_POSITION)
        {
            junctionReadPosition = junctionPositionFromMate(placement, up, down, upMatchedRead, primaryIsUpstream);
        }

        if(junctionReadPosition == NO_JUNCTION_POSITION)
        {
            if(mConfig.AnnotatedOnly)
            {
                return MergeOutcome.reject(RejectReason.NOVEL_JUNCTION);
            }

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
        ChrBaseRegion intron = intronAt(placement, up, down, upMatchedRead, junctionReadPosition);
        return MergeOutcome.success(up.Start, merged, intron, supp, primaryIsUpstream, mSpliceJunctions.strand(intron));
    }

    private static ChrBaseRegion intronAt(
            final Placement placement, final Side up, final Side down,
            final int upMatchedRead, final int readPosition)
    {
        int upLoss = upMatchedRead - readPosition;
        int downLoss = readPosition - down.LeadingS;
        return new ChrBaseRegion(
                placement.chromosome(), up.RefEnd - upLoss + 1, down.Start + downLoss - 1);
    }

    // Choose the highest junction tier; break ties deterministically per read.
    private int scanJunctionPositions(
            final Placement placement, final Side up, final Side down,
            final int upMatchedRead, final boolean primaryIsUpstream)
    {
        SpliceJunctions.Tier bestTier = SpliceJunctions.Tier.NONE;
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

            SpliceJunctions.Tier tier = mSpliceJunctions.tier(intronAt(placement, up, down, upMatchedRead, readPosition));
            if(tier == SpliceJunctions.Tier.NONE)
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

        int index = bestPositions.size() == 1 ? 0 : Math.floorMod(tieBreakSeed(placement, up), bestPositions.size());
        return bestPositions.get(index);
    }

    private static int tieBreakSeed(final Placement placement, final Side up)
    {
        int base = placement.readBases() != null ? Arrays.hashCode(placement.readBases()) : 0;
        return 31 * base + up.RefEnd;
    }

    // A mate hint pins the corresponding intron boundary.
    private int junctionPositionFromMate(
            final Placement placement, final Side up, final Side down,
            final int upMatchedRead, final boolean primaryIsUpstream)
    {
        for(ChrBaseRegion hint : placement.mateHintIntrons())
        {
            if(!hint.Chromosome.equals(placement.chromosome()))
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
        boolean RightExtend;
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
