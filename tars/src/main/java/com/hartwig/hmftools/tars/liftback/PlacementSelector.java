package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.XS_ATTRIBUTE;
import static com.hartwig.hmftools.tars.common.TarsConstants.CONFIDENT_MAPQ;
import static com.hartwig.hmftools.tars.common.TarsConstants.LOCAL_SV_MAX_LENGTH;
import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

import com.hartwig.hmftools.tars.common.ContigEntry;
import com.hartwig.hmftools.tars.liftback.features.OverhangGate;

import htsjdk.samtools.SAMRecord;

// Lifts a SAMRecord's alignments to genomic coordinates and decides the primary alignment, locus count and MAPQ.
// Every input record produces exactly one result.
public class PlacementSelector
{
    private final AlignmentLifter mAlignmentLifter;

    // when present, resolves hidden ties (XS==AS, no XA) on ref-only primaries landing inside an annotated exon.
    private final EnsemblAnnotationIndex mEnsemblAnnotationIndex;

    public PlacementSelector(final List<ContigEntry> entries)
    {
        this(entries, null);
    }

    public PlacementSelector(final List<ContigEntry> entries, final EnsemblAnnotationIndex annotationIndex)
    {
        mAlignmentLifter = new AlignmentLifter(entries);
        mEnsemblAnnotationIndex = annotationIndex;
    }

    ContigTranslator contigTranslator()
    {
        return mAlignmentLifter.translator();
    }

    // Convenience for non-discriminating callers: supplementaries, unmapped, lift-only paths, tests.
    public LiftedRecord resolve(final SAMRecord record)
    {
        return resolve(record, null, null);
    }

    public LiftedRecord resolve(final SAMRecord record, final OverhangGate overhangGate, final LiftedRecord mate)
    {
        if(record.getReadUnmappedFlag())
        {
            return LiftedRecord.unmapped("");
        }

        if(record.getSupplementaryAlignmentFlag())
        {
            return liftSupplementaryAlignment(record, overhangGate);
        }

        LiftedRecord alignments = liftPrimaryAlignments(record, overhangGate);
        if(!alignments.hasPlacement())
        {
            return alignments;
        }
        return selectPrimaryAlignment(record, alignments.liftedAlignments(), mate);
    }

    // README Steps 0-1: lift the primary and XA alts to genomic coordinates, then apply the overhang gate.
    public LiftedRecord liftPrimaryAlignments(final SAMRecord record, final OverhangGate overhangGate)
    {
        return mAlignmentLifter.liftPrimary(record, overhangGate);
    }

    public LiftedRecord selectPrimaryAlignment(
            final SAMRecord record, final List<LiftedAlignment> allAlignments, final LiftedRecord mate)
    {
        Selection selection = choosePrimaryAlignment(record, allAlignments, mate);
        return selectPrimaryAlignment(record, allAlignments, selection);
    }

    Selection choosePrimaryAlignment(
            final SAMRecord record, final List<LiftedAlignment> allAlignments, final LiftedRecord mate)
    {
        boolean hasMergeableSupplementary = hasMergeableSupplementary(allAlignments);
        return select(
                allAlignments, isConcordant(allAlignments), allAlignments.get(0), readSeed(record.getReadName()),
                record.getMappingQuality() != 0 && !hasMergeableSupplementary, mate);
    }

    PairSelection chooseMatePair(
            final SAMRecord first, final List<LiftedAlignment> firstAlignments,
            final SAMRecord second, final List<LiftedAlignment> secondAlignments)
    {
        boolean firstBwaPriority = first.getMappingQuality() != 0 && !hasMergeableSupplementary(firstAlignments);
        boolean secondBwaPriority = second.getMappingQuality() != 0 && !hasMergeableSupplementary(secondAlignments);
        return selectPair(
                firstAlignments, isConcordant(firstAlignments), firstAlignments.get(0), firstBwaPriority,
                first.getMappingQuality() == 0,
                secondAlignments, isConcordant(secondAlignments), secondAlignments.get(0), secondBwaPriority,
                second.getMappingQuality() == 0,
                readSeed(first.getReadName()));
    }

    LiftedRecord selectPrimaryAlignment(
            final SAMRecord record, final List<LiftedAlignment> allAlignments, final Selection selection)
    {
        int inputMapQuality = record.getMappingQuality();
        LiftedAlignment effectivePrimary = selection.alignment();

        for(LiftedAlignment alignment : allAlignments)
        {
            if(alignment != effectivePrimary && alignment.hasSupplementaryMerge())
            {
                alignment.Dropped = true;
            }
        }

        List<LiftedAlignment> keptAlignments = new ArrayList<>(allAlignments.size());
        for(LiftedAlignment alignment : allAlignments)
        {
            if(!alignment.Dropped)
            {
                keptAlignments.add(alignment);
            }
        }

        // Single kept alignment: the locus scan collapses to 1, so skip the per-read map.
        int numLoci = keptAlignments.size() == 1 ? 1 : countDistinctLoci(keptAlignments, effectivePrimary);

        boolean hiddenTie = inputMapQuality == 0 && hasHiddenTie(record);
        boolean inAnnotatedExon = mEnsemblAnnotationIndex != null
                && mEnsemblAnnotationIndex.containsExon(effectivePrimary.LiftedChromosome, effectivePrimary.LiftedPos);
        boolean randomTie = selection.reason().equals("random");
        int updatedMapQuality = decidePrimaryMapQuality(
                inputMapQuality, numLoci, hiddenTie, effectivePrimary.FromTxContig, inAnnotatedExon, randomTie);
        String note = selection.reason();

        if(effectivePrimary.hasSupplementaryMerge())
        {
            updatedMapQuality = Math.max(updatedMapQuality, effectivePrimary.MergedSupplementaryMapQuality);
            if(numLoci == 1)
            {
                updatedMapQuality = CONFIDENT_MAPQ;
            }
            note = appendNote(note, "supplementary-resolved");
        }

        if(selection.alignmentIndex() != 0)
        {
            TARS_LOGGER.trace(
                    "placement selection {}: primary -> {}:{} {} ({})",
                    record.getReadName(), effectivePrimary.LiftedChromosome, effectivePrimary.LiftedPos,
                    effectivePrimary.LiftedCigar, selection.reason());
        }

        return new LiftedRecord(updatedMapQuality, numLoci, note, selection.alignmentIndex(), allAlignments);
    }

    public static boolean usesGenomicScore(final List<LiftedAlignment> alignments, final int inputMapQuality)
    {
        return !isConcordant(alignments) && (inputMapQuality == 0 || hasMergeableSupplementary(alignments));
    }

    private static boolean hasMergeableSupplementary(final List<LiftedAlignment> alignments)
    {
        for(LiftedAlignment alignment : alignments)
        {
            if(alignment.hasSupplementaryMerge())
            {
                return true;
            }
        }
        return false;
    }

    private static String appendNote(final String existing, final String note)
    {
        if(existing == null || existing.isEmpty())
        {
            return note;
        }
        return existing + ";" + note;
    }

    // Lift a supplementary's own placement and XA alternatives together. The placement selected relative to the final
    // primary is used both for merging and for emission when the supplementary is not absorbed.
    public LiftedRecord liftSupplementaryAlignment(final SAMRecord record, final OverhangGate overhangGate)
    {
        return mAlignmentLifter.liftSupplementary(record, overhangGate);
    }

    // Distinct genomic loci among kept alignments: an alt overlapping the primary collapses into it; non-overlapping
    // alts interval-merge among themselves but are never chained back through the primary.
    private static int countDistinctLoci(final List<LiftedAlignment> alignments, final LiftedAlignment primary)
    {
        Map<String, List<int[]>> distinctSpans = new HashMap<>();
        distinctSpans.computeIfAbsent(primary.LiftedChromosome, k -> new ArrayList<>())
                .add(new int[] { primary.LiftedPos, primary.alignedEnd() });

        for(LiftedAlignment alignment : alignments)
        {
            if(alignment == primary || alignment.overlaps(primary))
            {
                continue;
            }
            distinctSpans.computeIfAbsent(alignment.LiftedChromosome, k -> new ArrayList<>())
                    .add(new int[] { alignment.LiftedPos, alignment.alignedEnd() });
        }

        int loci = 0;
        for(List<int[]> spans : distinctSpans.values())
        {
            spans.sort(Comparator.comparingInt(s -> s[0]));
            int clusterEnd = -1;
            for(int[] span : spans)
            {
                if(span[0] > clusterEnd)
                {
                    ++loci;
                    clusterEnd = span[1];
                }
                else
                {
                    clusterEnd = Math.max(clusterEnd, span[1]);
                }
            }
        }
        return loci;
    }

    // Emit-time NH recompute: drops Dropped alts and collapses alts overlapping the primary, so NH stays consistent
    // with the final XA. A placement-less record maps to one locus by definition.
    public static int countDistinctLoci(final LiftedRecord liftedRecord)
    {
        if(!liftedRecord.hasPlacement())
        {
            return 1;
        }

        List<LiftedAlignment> kept = new ArrayList<>(liftedRecord.liftedAlignments().size());
        for(LiftedAlignment alignment : liftedRecord.liftedAlignments())
        {
            if(!alignment.Dropped)
            {
                kept.add(alignment);
            }
        }
        return Math.max(countDistinctLoci(kept, liftedRecord.primaryAlignment()), 1);
    }

    // Only a MAPQ-0 primary is ever bumped, and only when it lifts to a single locus off a decisive pick. A hidden tie
    // (XS==AS, an equal-scoring alt bwa did not emit) blocks the bump unless tx provenance or an annotated exon vouches
    // for the placement. Anything bwa graded is left alone.
    static int decidePrimaryMapQuality(
            final int inputMapQuality, final int numLoci, final boolean hiddenTie,
            final boolean primaryFromTxContig, final boolean primaryInAnnotatedExon, final boolean randomTie)
    {
        boolean confident = inputMapQuality == 0 && !randomTie && numLoci == 1
                && (!hiddenTie || primaryFromTxContig || primaryInAnnotatedExon);

        return confident ? CONFIDENT_MAPQ : inputMapQuality;
    }

    // Deterministic per-read seed for the random placement picks; mates share a read name, so a pair is placed together.
    static int readSeed(final String readName)
    {
        int hash = readName.hashCode();
        return hash ^ (hash >>> 16);
    }

    // When XS == AS, an equally-scoring alt was not emitted by bwa; flag as a hidden tie to skip the MAPQ bump.
    private static boolean hasHiddenTie(final SAMRecord record)
    {
        Integer alignmentScore = record.getIntegerAttribute(ALIGNMENT_SCORE_ATTRIBUTE);
        Integer suboptimalScore = record.getIntegerAttribute(XS_ATTRIBUTE);
        return alignmentScore != null && suboptimalScore != null && suboptimalScore.intValue() == alignmentScore.intValue();
    }

    // Ref and tx agree on one contiguous placement, so there is nothing to choose between and the pick keeps bwa's
    // primary. An alt the overhang gate collapsed to a contiguous alignment is marked Dropped before placement selection
    // runs; it is a fabricated placement, so it contributes neither a source nor a locus.
    public static boolean isConcordant(final List<LiftedAlignment> alignments)
    {
        Set<AlignmentKey.Locus> loci = new HashSet<>();
        Set<String> distinctCigars = new HashSet<>();
        boolean hasRef = false;
        boolean hasTx = false;

        for(LiftedAlignment alignment : alignments)
        {
            if(alignment.Dropped)
            {
                continue;
            }

            // a surviving N means the two views disagree about splicing
            if(alignment.cigarHasN())
            {
                return false;
            }

            loci.add(alignment.key().locus());
            distinctCigars.add(alignment.LiftedCigar);

            if(alignment.FromTxContig)
            {
                hasTx = true;
            }
            else
            {
                hasRef = true;
            }
        }

        return hasRef && hasTx && loci.size() == 1 && distinctCigars.size() == 1;
    }

    // alignmentIndex is the winner's position in the alignment list, as stored by LiftedRecord.
    record Selection(int alignmentIndex, LiftedAlignment alignment, String reason)
    {
    }

    record PairSelection(Selection first, Selection second)
    {
    }

    // MAPQ-0 mates are selected as one fragment. Prefer the shortest DEL, then DUP, then INV within 1 Mb; pairs outside
    // that window are equivalent. Both placements are returned together, so read order cannot influence the result.
    static PairSelection selectPair(
            final List<LiftedAlignment> firstAlignments, final boolean firstConcordant,
            final LiftedAlignment firstSelf, final boolean firstBwaHasPriority,
            final List<LiftedAlignment> secondAlignments, final boolean secondConcordant,
            final LiftedAlignment secondSelf, final boolean secondBwaHasPriority, final int seed)
    {
        return selectPair(
                firstAlignments, firstConcordant, firstSelf, firstBwaHasPriority, !firstBwaHasPriority,
                secondAlignments, secondConcordant, secondSelf, secondBwaHasPriority, !secondBwaHasPriority, seed);
    }

    static PairSelection selectPair(
            final List<LiftedAlignment> firstAlignments, final boolean firstConcordant,
            final LiftedAlignment firstSelf, final boolean firstBwaHasPriority, final boolean firstMapqZero,
            final List<LiftedAlignment> secondAlignments, final boolean secondConcordant,
            final LiftedAlignment secondSelf, final boolean secondBwaHasPriority, final boolean secondMapqZero,
            final int seed)
    {
        Selection independentFirst = select(
                firstAlignments, firstConcordant, firstSelf, seed, firstBwaHasPriority);
        Selection independentSecond = select(
                secondAlignments, secondConcordant, secondSelf, seed, secondBwaHasPriority);

        Optional<DiscordantPairSelector.Choice> selection = DiscordantPairSelector.select(
                new DiscordantPairSelector.MatePlacements(
                        firstAlignments, independentFirst.alignment(), firstMapqZero && !firstConcordant),
                new DiscordantPairSelector.MatePlacements(
                        secondAlignments, independentSecond.alignment(), secondMapqZero && !secondConcordant),
                seed);
        if(selection.isEmpty())
        {
            return new PairSelection(independentFirst, independentSecond);
        }
        DiscordantPairSelector.Choice winner = selection.get();
        return new PairSelection(
                new Selection(indexOf(firstAlignments, winner.first()), winner.first(), winner.note()),
                new Selection(indexOf(secondAlignments, winner.second()), winner.second(), winner.note()));
    }

    // Mate-agnostic overload: single-end reads and callers with no lifted mate.
    static Selection select(
            final List<LiftedAlignment> alignments, final boolean concordant, final LiftedAlignment self,
            final int seed, final boolean bwaHasPriority)
    {
        return select(alignments, concordant, self, seed, bwaHasPriority, null);
    }

    // Returns the winning placement, its index and a short note. With bwaHasPriority false (MAPQ 0) placements are ranked by
    // recomputed genomic score, falling back to mate-proximity / junction / seed tie-breaks only on a score tie or unscored
    // placements (split read left for supplementary-resolve). bwaHasPriority true leaves bwa's order untouched.
    static Selection select(
            final List<LiftedAlignment> alignments, final boolean concordant, final LiftedAlignment self,
            final int seed, final boolean bwaHasPriority, final LiftedRecord mate)
    {
        if(concordant || bwaHasPriority)
        {
            return keepBwaPrimary(alignments, self);
        }
        return pickByScore(alignments, self, seed, mate);
    }

    // Highest recomputed genome score wins ("score"). Top-score ties are settled in order by: closest plausible mate,
    // supplementary support, junction over soft clip, then a read-name seed. Nothing is dropped; losers ride in XA.
    private static Selection pickByScore(
            final List<LiftedAlignment> alignments, final LiftedAlignment self, final int seed, final LiftedRecord mate)
    {
        List<LiftedAlignment> placements = new ArrayList<>();
        for(LiftedAlignment alignment : alignments)
        {
            if(!alignment.Dropped)
            {
                placements.add(alignment);
            }
        }
        if(placements.size() < 2)
        {
            return keepBwaPrimary(alignments, self);
        }

        int topScore = Integer.MIN_VALUE;
        for(LiftedAlignment alignment : placements)
        {
            topScore = Math.max(topScore, alignment.GenomicScore);
        }
        if(topScore == Integer.MIN_VALUE && !hasMergeableSupplementary(placements))
        {
            return keepBwaPrimary(alignments, self);
        }

        // Collapse identical placements (same locus + CIGAR from different sources, e.g. a ref self and a tx alt that
        // lift to the same contiguous alignment) so the tie is over distinct placements, not weighted by source count.
        List<LiftedAlignment> top = new ArrayList<>();
        Set<AlignmentKey> topKeys = new HashSet<>();
        for(LiftedAlignment alignment : placements)
        {
            if(alignment.GenomicScore == topScore && topKeys.add(alignment.key()))
            {
                top.add(alignment);
            }
        }

        boolean tie = top.size() > 1;
        LiftedAlignment winner;
        String note;
        if(!tie)
        {
            winner = top.get(0);
            note = "score";
        }
        else
        {
            List<LiftedAlignment> contenders = closestMateSubset(top, mate);
            if(contenders.size() == 1)
            {
                winner = contenders.get(0);
                note = "mate";
            }
            else
            {
                contenders = supplementarySupportedSubset(contenders);
                if(contenders.size() == 1)
                {
                    winner = contenders.get(0);
                    note = "supplementary";
                }
                else
                {
                    LiftedAlignment junction = preferJunctionOverSoftClip(contenders);
                    if(junction != null)
                    {
                        winner = junction;
                        note = "junction";
                    }
                    else
                    {
                        winner = contenders.get(Math.floorMod(seed, contenders.size()));
                        note = "random";
                    }
                }
            }
        }
        return new Selection(indexOf(alignments, winner), winner, note);
    }

    private static List<LiftedAlignment> supplementarySupportedSubset(final List<LiftedAlignment> top)
    {
        List<LiftedAlignment> supported = new ArrayList<>();
        for(LiftedAlignment alignment : top)
        {
            if(alignment.hasSupplementaryMerge())
            {
                supported.add(alignment);
            }
        }
        return supported.isEmpty() ? top : supported;
    }

    private static Selection keepBwaPrimary(final List<LiftedAlignment> alignments, final LiftedAlignment self)
    {
        return new Selection(indexOf(alignments, self), self, "");
    }

    // identity, not equals: a placement revised by withLiftedCigar is a distinct object at the same list position.
    private static int indexOf(final List<LiftedAlignment> alignments, final LiftedAlignment target)
    {
        for(int i = 0; i < alignments.size(); ++i)
        {
            if(alignments.get(i) == target)
            {
                return i;
            }
        }
        return LiftedRecord.NO_PRIMARY;
    }

    // Prefer the candidate with the smallest gap to any viable mate placement. The 1 Mb limit prevents a distant
    // same-chromosome mate from influencing the decision; an absent mate, an out-of-range mate or a distance tie leaves
    // the full contender set unchanged.
    private static List<LiftedAlignment> closestMateSubset(final List<LiftedAlignment> top, final LiftedRecord mate)
    {
        if(mate == null || !mate.hasPlacement())
        {
            return top;
        }

        int closestDistance = Integer.MAX_VALUE;
        for(LiftedAlignment alignment : top)
        {
            closestDistance = Math.min(closestDistance, mateDistance(alignment, mate));
        }
        if(closestDistance > LOCAL_SV_MAX_LENGTH)
        {
            return top;
        }

        List<LiftedAlignment> closest = new ArrayList<>();
        for(LiftedAlignment alignment : top)
        {
            if(mateDistance(alignment, mate) == closestDistance)
            {
                closest.add(alignment);
            }
        }
        return closest.size() == top.size() ? top : closest;
    }

    private static int mateDistance(final LiftedAlignment alignment, final LiftedRecord mate)
    {
        int closestDistance = Integer.MAX_VALUE;
        for(LiftedAlignment mateAlignment : mate.liftedAlignments())
        {
            if(!mateAlignment.Dropped)
            {
                closestDistance = Math.min(closestDistance, mateDistance(alignment, mateAlignment));
            }
        }
        return closestDistance;
    }

    private static int mateDistance(final LiftedAlignment alignment, final LiftedAlignment mate)
    {
        return alignment.alignedBlockDistance(mate);
    }

    // Tie-break within an equal-top-score set: a spliced placement (real N junction) beats a clipped placement
    // (soft-clip, no N) at the same lifted locus. bwa soft-clipped rather than cross the intron, so the junction is
    // the correct RNA interpretation and is not left to the seed.
    private static LiftedAlignment preferJunctionOverSoftClip(final List<LiftedAlignment> top)
    {
        for(LiftedAlignment junction : top)
        {
            if(!junction.cigarHasN())
            {
                continue;
            }
            for(LiftedAlignment clipped : top)
            {
                if(clipped != junction
                        && clipped.locusKey().equals(junction.locusKey())
                        && clipped.cigarHasSoftClip()
                        && !clipped.cigarHasN())
                {
                    return junction;
                }
            }
        }
        return null;
    }

}
