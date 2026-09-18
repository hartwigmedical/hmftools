package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.XS_ATTRIBUTE;
import static com.hartwig.hmftools.tars.common.TarsConstants.CONFIDENT_MAPQ;
import static com.hartwig.hmftools.tars.common.TarsConstants.LOCAL_SV_MAX_LENGTH;
import static com.hartwig.hmftools.tars.common.TarsConstants.PRIMARY_AS_UNMAP_THRESHOLD;
import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.tars.common.ContigEntry;
import com.hartwig.hmftools.tars.liftback.features.OverhangGate;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public final class AlignmentSelector
{
    private static final long LOCAL_PAIR_MAX_SEPARATION = 1000;
    private static final int LOCAL = 0;
    private static final int DELETION = 1;
    private static final int DUPLICATION = 2;
    private static final int INVERSION = 3;
    private static final int OTHER = 4;
    private static final PairPriority FALLBACK_PRIORITY = new PairPriority(OTHER, Long.MAX_VALUE);

    private final AlignmentLifter mAlignmentLifter;

    private final EnsemblAnnotationIndex mEnsemblAnnotationIndex;

    AlignmentSelector(final List<ContigEntry> entries)
    {
        this(entries, null);
    }

    AlignmentSelector(final List<ContigEntry> entries, final EnsemblAnnotationIndex annotationIndex)
    {
        mAlignmentLifter = new AlignmentLifter(entries);
        mEnsemblAnnotationIndex = annotationIndex;
    }

    ContigTranslator contigTranslator()
    {
        return mAlignmentLifter.translator();
    }

    LiftedRecord resolve(final SAMRecord record)
    {
        return resolve(record, null, null);
    }

    LiftedRecord resolve(final SAMRecord record, final OverhangGate overhangGate, final LiftedRecord mate)
    {
        if(record.getReadUnmappedFlag())
        {
            return LiftedRecord.unmapped("");
        }

        if(record.getSupplementaryAlignmentFlag())
        {
            return finaliseSupplementary(liftSupplementaryAlignment(record, overhangGate));
        }

        LiftedRecord alignments = liftPrimaryAlignments(record, overhangGate);
        if(!alignments.hasPlacement())
        {
            return alignments;
        }
        return selectPrimaryAlignment(record, alignments.liftedAlignments(), mate);
    }

    LiftedRecord liftPrimaryAlignments(final SAMRecord record, final OverhangGate overhangGate)
    {
        return mAlignmentLifter.liftPrimary(record, overhangGate);
    }

    LiftedRecord selectPrimaryAlignment(
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

        dropUnselectedMerges(allAlignments, effectivePrimary);
        List<LiftedAlignment> keptAlignments = collectKeptAlignments(allAlignments);
        traceCandidates(record, keptAlignments, effectivePrimary);

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

        traceScoreRegression(record, keptAlignments, effectivePrimary, selection.reason());

        return new LiftedRecord(updatedMapQuality, numLoci, note, selection.alignmentIndex(), allAlignments);
    }

    private static void dropUnselectedMerges(
            final List<LiftedAlignment> alignments, final LiftedAlignment selected)
    {
        for(LiftedAlignment alignment : alignments)
        {
            if(alignment != selected && alignment.hasSupplementaryMerge())
            {
                alignment.Dropped = true;
            }
        }
    }

    private static List<LiftedAlignment> collectKeptAlignments(final List<LiftedAlignment> alignments)
    {
        List<LiftedAlignment> kept = new ArrayList<>(alignments.size());
        for(LiftedAlignment alignment : alignments)
        {
            if(!alignment.Dropped)
            {
                kept.add(alignment);
            }
        }
        return kept;
    }

    private static void traceCandidates(
            final SAMRecord record, final List<LiftedAlignment> alignments, final LiftedAlignment selected)
    {
        if(!TARS_LOGGER.isTraceEnabled() || alignments.size() < 2)
        {
            return;
        }

        for(LiftedAlignment alignment : alignments)
        {
            TARS_LOGGER.trace(
                    "candidate {}: mapQuality={} {}:{} {} score={} fromTx={} chosen={}",
                    record.getReadName(), record.getMappingQuality(),
                    alignment.LiftedChromosome, alignment.LiftedPos, alignment.LiftedCigar,
                    alignment.GenomicScore, alignment.FromTxContig, alignment == selected);
        }
    }

    private static void traceScoreRegression(
            final SAMRecord record, final List<LiftedAlignment> alignments,
            final LiftedAlignment selected, final String reason)
    {
        if(!TARS_LOGGER.isTraceEnabled())
        {
            return;
        }

        LiftedAlignment bestScored = selected;
        for(LiftedAlignment alignment : alignments)
        {
            if(alignment.GenomicScore > bestScored.GenomicScore)
            {
                bestScored = alignment;
            }
        }

        if(bestScored != selected)
        {
            TARS_LOGGER.trace(
                    "score regression {}: inputMapQuality={} reason={} selected {}:{} {} aligned={} score={}"
                            + " over {}:{} {} aligned={} score={}",
                    record.getReadName(), record.getMappingQuality(), reason,
                    selected.LiftedChromosome, selected.LiftedPos, selected.LiftedCigar,
                    alignedBaseCount(selected.LiftedCigar), selected.GenomicScore,
                    bestScored.LiftedChromosome, bestScored.LiftedPos, bestScored.LiftedCigar,
                    alignedBaseCount(bestScored.LiftedCigar), bestScored.GenomicScore);
        }
    }

    private static int alignedBaseCount(final String cigar)
    {
        if(cigar == null)
        {
            return 0;
        }

        int alignedBases = 0;
        for(CigarElement element : CigarUtils.cigarFromStr(cigar).getCigarElements())
        {
            if(element.getOperator().isAlignment())
            {
                alignedBases += element.getLength();
            }
        }
        return alignedBases;
    }

    static boolean usesGenomicScore(final List<LiftedAlignment> alignments, final int inputMapQuality)
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

    LiftedRecord liftSupplementaryAlignment(final SAMRecord record, final OverhangGate overhangGate)
    {
        return mAlignmentLifter.liftSupplementary(record, overhangGate);
    }

    static LiftedRecord finaliseSupplementary(final LiftedRecord alignment)
    {
        if(alignment.hasPlacement() && alignment.updatedMapQuality() == 0 && alignment.xaTag() == null)
        {
            return alignment.withMapQuality(CONFIDENT_MAPQ);
        }
        return alignment;
    }

    // Alternatives collapse into the primary or each other, but do not chain through the primary.
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

    static int countDistinctLoci(final LiftedRecord liftedRecord)
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

    // XS == AS blocks a MAPQ bump unless transcript evidence resolves the hidden tie.
    static int decidePrimaryMapQuality(
            final int inputMapQuality, final int numLoci, final boolean hiddenTie,
            final boolean primaryFromTxContig, final boolean primaryInAnnotatedExon, final boolean randomTie)
    {
        boolean confident = inputMapQuality == 0 && !randomTie && numLoci == 1
                && (!hiddenTie || primaryFromTxContig || primaryInAnnotatedExon);

        return confident ? CONFIDENT_MAPQ : inputMapQuality;
    }

    static int readSeed(final String readName)
    {
        int hash = readName.hashCode();
        return hash ^ (hash >>> 16);
    }

    private static boolean hasHiddenTie(final SAMRecord record)
    {
        Integer alignmentScore = record.getIntegerAttribute(ALIGNMENT_SCORE_ATTRIBUTE);
        Integer suboptimalScore = record.getIntegerAttribute(XS_ATTRIBUTE);
        return alignmentScore != null && suboptimalScore != null && suboptimalScore.intValue() == alignmentScore.intValue();
    }

    static boolean isConcordant(final List<LiftedAlignment> alignments)
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

    record Selection(int alignmentIndex, LiftedAlignment alignment, String reason)
    {
    }

    record PairSelection(Selection first, Selection second)
    {
    }

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

        Optional<PairChoice> selection = selectMateAlignments(
                new MateAlignments(
                        firstAlignments, independentFirst.alignment(), firstMapqZero && !firstConcordant),
                new MateAlignments(
                        secondAlignments, independentSecond.alignment(), secondMapqZero && !secondConcordant),
                seed);
        if(selection.isEmpty())
        {
            return new PairSelection(independentFirst, independentSecond);
        }
        PairChoice winner = selection.get();
        return new PairSelection(
                new Selection(indexOf(firstAlignments, winner.first()), winner.first(), winner.note()),
                new Selection(indexOf(secondAlignments, winner.second()), winner.second(), winner.note()));
    }

    private static Optional<PairChoice> selectMateAlignments(
            final MateAlignments first, final MateAlignments second, final int seed)
    {
        List<AlignmentPair> bestPairs = new ArrayList<>();
        PairPriority bestPriority = null;

        for(LiftedAlignment firstAlignment : mateCandidates(first))
        {
            for(LiftedAlignment secondAlignment : mateCandidates(second))
            {
                PairPriority priority = matePairPriority(firstAlignment, secondAlignment);
                AlignmentPair pair = new AlignmentPair(firstAlignment, secondAlignment);
                int comparison = bestPriority == null ? -1 : priority.compareTo(bestPriority);
                if(comparison < 0)
                {
                    bestPairs.clear();
                    bestPriority = priority;
                }
                if(comparison <= 0)
                {
                    bestPairs.add(pair);
                }
            }
        }

        if(bestPairs.isEmpty())
        {
            return Optional.empty();
        }

        List<AlignmentPair> contenders = topScoringPairs(bestPairs);
        contenders.sort(Comparator.comparing(AlignmentPair::canonicalKey));
        AlignmentPair winner = contenders.get(Math.floorMod(seed, contenders.size()));
        String note = bestPairs.size() == 1 ? "mate" : (contenders.size() == 1 ? "score" : "random");
        return Optional.of(new PairChoice(winner.first(), winner.second(), note));
    }

    private static List<LiftedAlignment> mateCandidates(final MateAlignments mate)
    {
        if(!mate.selectable())
        {
            return List.of(mate.selected());
        }

        List<LiftedAlignment> candidates = new ArrayList<>();
        Set<AlignmentKey> seen = new HashSet<>();
        for(LiftedAlignment alignment : mate.alignments())
        {
            if(!alignment.Dropped
                    && (alignment.GenomicScore == Integer.MIN_VALUE
                            || alignment.GenomicScore >= PRIMARY_AS_UNMAP_THRESHOLD)
                    && seen.add(alignment.key()))
            {
                candidates.add(alignment);
            }
        }
        if(!candidates.isEmpty())
        {
            return candidates;
        }

        for(LiftedAlignment alignment : mate.alignments())
        {
            if(!alignment.Dropped && seen.add(alignment.key()))
            {
                candidates.add(alignment);
            }
        }
        return candidates;
    }

    private static List<AlignmentPair> topScoringPairs(final List<AlignmentPair> pairs)
    {
        AlignmentPair firstPair = pairs.get(0);
        boolean fixedFirst = true;
        boolean fixedSecond = true;
        for(AlignmentPair pair : pairs)
        {
            fixedFirst &= pair.first() == firstPair.first();
            fixedSecond &= pair.second() == firstPair.second();
        }

        long topScore = Long.MIN_VALUE;
        for(AlignmentPair pair : pairs)
        {
            if((!fixedFirst && pair.first().GenomicScore == Integer.MIN_VALUE)
                    || (!fixedSecond && pair.second().GenomicScore == Integer.MIN_VALUE))
            {
                return pairs;
            }
            topScore = Math.max(topScore, pairScore(pair, fixedFirst, fixedSecond));
        }

        List<AlignmentPair> topScored = new ArrayList<>();
        for(AlignmentPair pair : pairs)
        {
            if(pairScore(pair, fixedFirst, fixedSecond) == topScore)
            {
                topScored.add(pair);
            }
        }
        return topScored;
    }

    private static long pairScore(
            final AlignmentPair pair, final boolean fixedFirst, final boolean fixedSecond)
    {
        long score = fixedFirst ? 0 : pair.first().GenomicScore;
        return score + (fixedSecond ? 0 : pair.second().GenomicScore);
    }

    private static PairPriority matePairPriority(
            final LiftedAlignment first, final LiftedAlignment second)
    {
        int firstBreakend = first.ForwardStrand ? first.alignedEnd() : first.LiftedPos;
        int secondBreakend = second.ForwardStrand ? second.alignedEnd() : second.LiftedPos;
        Orientation firstOrientation = first.ForwardStrand ? Orientation.FORWARD : Orientation.REVERSE;
        Orientation secondOrientation = second.ForwardStrand ? Orientation.FORWARD : Orientation.REVERSE;

        long separation = first.alignedBlockDistance(second);
        if(first.LiftedChromosome.equals(second.LiftedChromosome)
                && separation <= LOCAL_PAIR_MAX_SEPARATION)
        {
            return new PairPriority(LOCAL, separation);
        }

        return pairPriority(
                first.LiftedChromosome, firstBreakend, firstOrientation,
                second.LiftedChromosome, secondBreakend, secondOrientation,
                separation);
    }

    private record MateAlignments(
            List<LiftedAlignment> alignments, LiftedAlignment selected, boolean selectable)
    {
        private MateAlignments
        {
            alignments = List.copyOf(alignments);
        }
    }

    private record PairChoice(LiftedAlignment first, LiftedAlignment second, String note) { }

    private record AlignmentPair(LiftedAlignment first, LiftedAlignment second)
    {
        String canonicalKey()
        {
            String firstKey = first.key().toString();
            String secondKey = second.key().toString();
            return firstKey.compareTo(secondKey) <= 0
                    ? firstKey + '|' + secondKey
                    : secondKey + '|' + firstKey;
        }
    }

    static List<RecordAlignment> selectSupplementaryAlignments(
            final String primaryChromosome, final boolean primaryForwardStrand,
            final int primaryStart, final String primaryCigar, final byte[] readBases,
            final List<RecordAlignment> supplementaries)
    {
        return selectSupplementaryAlignments(
                primaryChromosome, primaryForwardStrand, primaryStart,
                CigarUtils.cigarElementsFromStr(primaryCigar), readBases, supplementaries);
    }

    public static List<RecordAlignment> selectSupplementaryAlignments(
            final String primaryChromosome, final boolean primaryForwardStrand,
            final int primaryStart, final List<CigarElement> primaryCigar, final byte[] readBases,
            final List<RecordAlignment> supplementaries)
    {
        Map<Integer, RecordAlignment> selected = new LinkedHashMap<>();

        for(RecordAlignment candidate : supplementaries)
        {
            RecordAlignment current = selected.get(candidate.recordIndex());
            if(current == null)
            {
                selected.put(candidate.recordIndex(), candidate);
            }
            else if(current.mapQuality() == 0
                    && compareSupplementaryAlignments(
                            primaryChromosome, primaryForwardStrand, primaryStart, primaryCigar,
                            readBases, candidate, current) < 0)
            {
                selected.put(candidate.recordIndex(), candidate);
            }
        }

        return new ArrayList<>(selected.values());
    }

    private static int compareSupplementaryAlignments(
            final String primaryChromosome, final boolean primaryForwardStrand,
            final int primaryStart, final List<CigarElement> primaryCigar, final byte[] readBases,
            final RecordAlignment candidate, final RecordAlignment current)
    {
        PairPriority candidatePriority = supplementaryPriority(
                primaryChromosome, primaryForwardStrand, primaryStart, primaryCigar, candidate);
        PairPriority currentPriority = supplementaryPriority(
                primaryChromosome, primaryForwardStrand, primaryStart, primaryCigar, current);

        int priorityComparison = candidatePriority.compareTo(currentPriority);
        if(priorityComparison != 0)
        {
            return priorityComparison;
        }

        int candidateOrder = 31 * Arrays.hashCode(readBases) + supplementaryHash(candidate);
        int currentOrder = 31 * Arrays.hashCode(readBases) + supplementaryHash(current);
        return Integer.compareUnsigned(candidateOrder, currentOrder);
    }

    private static int supplementaryHash(final RecordAlignment supplementary)
    {
        LiftedAlignment alignment = supplementary.alignment();
        return new SupplementaryKey(
                supplementary.recordIndex(), alignment.LiftedChromosome, alignment.ForwardStrand,
                alignment.LiftedPos, alignment.LiftedCigar, supplementary.mapQuality()).hashCode();
    }

    private static PairPriority supplementaryPriority(
            final String primaryChromosome, final boolean primaryForwardStrand,
            final int primaryStart, final List<CigarElement> primaryCigar,
            final RecordAlignment supplementary)
    {
        LiftedAlignment alignment = supplementary.alignment();
        if(!primaryChromosome.equals(alignment.LiftedChromosome))
        {
            return FALLBACK_PRIORITY;
        }

        List<CigarElement> supplementaryCigar = CigarUtils.cigarElementsFromStr(alignment.LiftedCigar);
        int primaryLeadingClip = CigarUtils.leftSoftClipLength(primaryCigar);
        int supplementaryLeadingClip = CigarUtils.leftSoftClipLength(supplementaryCigar);
        if(primaryLeadingClip == supplementaryLeadingClip)
        {
            return FALLBACK_PRIORITY;
        }

        boolean primaryLinksEnd = primaryLeadingClip < supplementaryLeadingClip;
        boolean supplementaryLinksEnd = !primaryLinksEnd;
        int primaryEnd = primaryStart + CigarUtils.cigarAlignedLength(primaryCigar) - 1;
        int supplementaryEnd = alignment.LiftedPos + CigarUtils.cigarAlignedLength(supplementaryCigar) - 1;

        return pairPriority(
                primaryChromosome,
                breakendPosition(primaryStart, primaryEnd, primaryForwardStrand, primaryLinksEnd),
                breakendOrientation(primaryForwardStrand, primaryLinksEnd),
                alignment.LiftedChromosome,
                breakendPosition(
                        alignment.LiftedPos, supplementaryEnd,
                        alignment.ForwardStrand, supplementaryLinksEnd),
                breakendOrientation(alignment.ForwardStrand, supplementaryLinksEnd));
    }

    public record RecordAlignment(
            int recordIndex, int alignmentIndex, LiftedAlignment alignment, int mapQuality)
    {
    }

    private record SupplementaryKey(
            int index, String chromosome, boolean forwardStrand, int start, String cigar, int mapQuality)
    {
    }

    private static int breakendPosition(
            final int start, final int end, final boolean forwardStrand, final boolean linksEnd)
    {
        if(linksEnd)
        {
            return forwardStrand ? end : start;
        }
        return forwardStrand ? start : end;
    }

    private static Orientation breakendOrientation(final boolean forwardStrand, final boolean linksEnd)
    {
        return linksEnd == forwardStrand ? Orientation.FORWARD : Orientation.REVERSE;
    }

    private static PairPriority pairPriority(
            final String firstChromosome, final int firstPosition, final Orientation firstOrientation,
            final String secondChromosome, final int secondPosition, final Orientation secondOrientation)
    {
        return pairPriority(
                firstChromosome, firstPosition, firstOrientation,
                secondChromosome, secondPosition, secondOrientation,
                Math.abs((long) firstPosition - secondPosition));
    }

    private static PairPriority pairPriority(
            final String firstChromosome, final int firstPosition, final Orientation firstOrientation,
            final String secondChromosome, final int secondPosition, final Orientation secondOrientation,
            final long separation)
    {
        if(!firstChromosome.equals(secondChromosome) || separation > LOCAL_SV_MAX_LENGTH)
        {
            return FALLBACK_PRIORITY;
        }
        if(firstOrientation == secondOrientation)
        {
            return new PairPriority(INVERSION, separation);
        }

        long breakendLength = Math.abs((long) firstPosition - secondPosition);
        if(breakendLength == 0)
        {
            return new PairPriority(DUPLICATION, separation);
        }

        boolean firstIsLower = firstPosition < secondPosition;
        int category = firstIsLower == firstOrientation.isForward() ? DELETION : DUPLICATION;
        return new PairPriority(category, separation);
    }

    private record PairPriority(int category, long length) implements Comparable<PairPriority>
    {
        @Override
        public int compareTo(final PairPriority other)
        {
            int categoryComparison = Integer.compare(category, other.category);
            return categoryComparison != 0 ? categoryComparison : Long.compare(length, other.length);
        }
    }

    static Selection select(
            final List<LiftedAlignment> alignments, final boolean concordant, final LiftedAlignment self,
            final int seed, final boolean bwaHasPriority)
    {
        return select(alignments, concordant, self, seed, bwaHasPriority, null);
    }

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

        // Do not weight a placement by the number of sources that produced it.
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
