package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.XA_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.XS_ATTRIBUTE;
import static com.hartwig.hmftools.tars.common.TarsConstants.CONFIDENT_MAPQ;
import static com.hartwig.hmftools.tars.common.TarsConstants.MATE_PROXIMITY_MAX_DISTANCE;
import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.tars.common.ContigEntry;

import htsjdk.samtools.SAMRecord;

// Lifts a SAMRecord's alignments to genomic coordinates and decides its outcome: which candidate becomes the primary,
// how many loci it maps to, and what MAPQ that placement earns. Every input record produces exactly one result.
public class LiftBackDiscriminator
{
    private final ContigTranslator mContigTranslator;

    // when present, resolves hidden ties (XS==AS, no XA) on ref-only primaries landing inside an annotated exon.
    private final ExonRegionIndex mExonIndex;

    public LiftBackDiscriminator(final List<ContigEntry> entries)
    {
        this(entries, null);
    }

    public LiftBackDiscriminator(final List<ContigEntry> entries, final ExonRegionIndex exonIndex)
    {
        mContigTranslator = new ContigTranslator(entries);
        mExonIndex = exonIndex;
    }

    ContigTranslator contigTranslator()
    {
        return mContigTranslator;
    }

    private LiftedAlignment liftSelf(final SAMRecord record)
    {
        LiftedAlignment lifted = mContigTranslator.liftAlignment(
                record.getReferenceName(), record.getAlignmentStart(), record.getCigarString(),
                getInt(record, NUM_MUTATONS_ATTRIBUTE), !record.getReadNegativeStrandFlag());

        if(lifted == null)
        {
            logLiftFailure(record);
        }

        return lifted;
    }

    // A lift that produces nothing is rare, so log it per read against the transcript segment owning the position:
    // the cause reads off the numbers. pos below segStart or readEnd above segEnd is an overhang the clamp could not
    // absorb, pos above segEnd is an inter-transcript spacer hit, and wholly inside means the walk ran off the last exon.
    private void logLiftFailure(final SAMRecord record)
    {
        String contig = record.getReferenceName();
        int pos = record.getAlignmentStart();
        int readEnd = pos + record.getCigar().getReferenceLength() - 1;
        String role = record.getSupplementaryAlignmentFlag() ? "supp" : "primary";
        ContigEntry segment = mContigTranslator.findSegment(contig, pos);

        if(segment == null)
        {
            TARS_LOGGER.debug(
                    "lift failed {} {}: {}:{}-{} {} - contig has no segments",
                    record.getReadName(), role, contig, pos, readEnd, record.getCigarString());
            return;
        }

        TARS_LOGGER.debug(
                "lift failed {} {}: {}:{}-{} {} - segment {} [{}-{}] exons({})",
                record.getReadName(), role, contig, pos, readEnd, record.getCigarString(),
                segment.transName(), segment.contigStart(), segment.contigEnd(), segment.exonSpans().size());
    }

    // No-reconcile convenience for non-discriminating callers (supplementaries, unmapped, lift-only paths, tests).
    public LiftedRecord resolve(final SAMRecord record)
    {
        return resolve(record, null, null);
    }

    public LiftedRecord resolve(final SAMRecord record, final OverhangReconciler reconciler)
    {
        return resolve(record, reconciler, null);
    }

    public LiftedRecord resolve(final SAMRecord record, final OverhangReconciler reconciler, final LiftedRecord mate)
    {
        if(record.getReadUnmappedFlag())
        {
            return LiftedRecord.unmapped("");
        }

        if(record.getSupplementaryAlignmentFlag())
        {
            return liftSupplementary(record);
        }

        return resolvePrimary(record, reconciler, mate);
    }

    private LiftedRecord resolvePrimary(final SAMRecord record, final OverhangReconciler reconciler, final LiftedRecord mate)
    {
        LiftedAlignment self = liftSelf(record);

        if(self == null)
        {
            return LiftedRecord.unmapped("primary_translate_failed");
        }

        List<LiftedAlignment> alts = parseAndLiftXa(record);
        List<LiftedAlignment> allAlignments = new ArrayList<>(1 + alts.size());
        allAlignments.add(self);
        allAlignments.addAll(alts);

        // Normalize each candidate's cigar before discriminating so the ref-vs-tx features are measured, not assumed.
        // Null on the lift-only paths, where nothing is being decided.
        if(reconciler != null)
        {
            reconciler.reconcileAlignmentsToGenome(allAlignments, record);
        }

        // reconciling can replace self (index 0) with a revised copy, so re-fetch it before the pick.
        self = allAlignments.get(0);

        boolean concordant = isConcordant(allAlignments);
        int inputMapQuality = record.getMappingQuality();
        int seed = readSeed(record.getReadName());
        ApplyResult outcome = apply(
                allAlignments, concordant, self, seed, inputMapQuality != 0, mate);
        LiftedAlignment effectivePrimary = outcome.effectivePrimary();

        List<LiftedAlignment> keptAlignments = new ArrayList<>(allAlignments.size());
        for(LiftedAlignment alignment : allAlignments)
        {
            if(!alignment.Dropped)
            {
                keptAlignments.add(alignment);
            }
        }

        // Single kept candidate: the locus scan collapses to 1, so skip the per-read map.
        int numLoci = keptAlignments.size() == 1 ? 1 : countDistinctLoci(keptAlignments, effectivePrimary);

        boolean hiddenTie = inputMapQuality == 0 && hasHiddenTie(record);
        boolean inAnnotatedExon = mExonIndex != null
                && mExonIndex.contains(effectivePrimary.LiftedChromosome, effectivePrimary.LiftedPos);
        boolean randomTie = outcome.note().equals("random");
        int updatedMapQuality = decidePrimaryMapQuality(
                inputMapQuality, numLoci, hiddenTie, effectivePrimary.FromTxContig, inAnnotatedExon, randomTie);

        if(outcome.primaryIndex() != 0)
        {
            TARS_LOGGER.trace(
                    "discriminator {}: primary -> {}:{} {} ({})",
                    record.getReadName(), effectivePrimary.LiftedChromosome, effectivePrimary.LiftedPos,
                    effectivePrimary.LiftedCigar, outcome.note());
        }

        return new LiftedRecord(updatedMapQuality, numLoci, outcome.note(), outcome.primaryIndex(), allAlignments);
    }

    // A supplementary is only lifted, never discriminated: lift its own coords with no XA parse or locus pick.
    // MAPQ=0 on a tx-contig supplementary is the multi-alt-contig tie artefact; bump to CONFIDENT_MAPQ.
    private LiftedRecord liftSupplementary(final SAMRecord record)
    {
        LiftedAlignment lifted = liftSelf(record);

        if(lifted == null)
        {
            return LiftedRecord.unmapped("supp_translate_failed");
        }

        int inputMapQuality = record.getMappingQuality();
        int outputMapQuality = (lifted.FromTxContig && inputMapQuality == 0) ? CONFIDENT_MAPQ : inputMapQuality;

        return new LiftedRecord(outputMapQuality, 1, "", 0, List.of(lifted));
    }

    // Self is excluded from the dedup key set so a Tx XA alt lifting to the same coords as a ref self is preserved (drives CONCORDANT).
    private List<LiftedAlignment> parseAndLiftXa(final SAMRecord record)
    {
        List<LiftedAlignment> alts = new ArrayList<>();
        Set<String> seenKeys = new HashSet<>();

        for(LiftedAlignment lifted : mContigTranslator.liftXaAlignments(record.getStringAttribute(XA_ATTRIBUTE)))
        {
            if(seenKeys.add(liftedKey(lifted)))
            {
                alts.add(lifted);
            }
        }

        return alts;
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
    // with the XA after the soft-clip extension mutates the shared list. A placement-less record maps to one locus by definition.
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

    private static String liftedKey(final LiftedAlignment alignment)
    {
        return alignment.LiftedChromosome + ":" + alignment.LiftedPos + ":" + alignment.LiftedCigar
                + ":" + (alignment.ForwardStrand ? '+' : '-');
    }

    private static int getInt(final SAMRecord record, final String tag)
    {
        Integer value = record.getIntegerAttribute(tag);
        return value != null ? value : 0;
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
    // primary. An alt the overhang gate collapsed to a contiguous alignment is marked Dropped before the discriminator
    // runs; it is a fabricated placement, so it contributes neither a source nor a locus.
    public static boolean isConcordant(final List<LiftedAlignment> alignments)
    {
        Set<String> loci = new HashSet<>();
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

            loci.add(alignment.locusKey());
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

    // primaryIndex is the winner's position in the alignment list, which is what LiftedRecord stores.
    public record ApplyResult(int primaryIndex, LiftedAlignment effectivePrimary, String note)
    {
    }

    // Mate-agnostic overload: single-end reads and callers with no lifted mate.
    public static ApplyResult apply(
            final List<LiftedAlignment> alignments, final boolean concordant, final LiftedAlignment self,
            final int seed, final boolean bwaHasPriority)
    {
        return apply(alignments, concordant, self, seed, bwaHasPriority, null);
    }

    // Returns the winning candidate, its index and a short note. When bwa expressed no priority
    // (bwaHasPriority false, ie MAPQ 0) it ranks candidates by recomputed genomic score and takes the clear winner, falling
    // back to mate-proximity / junction / seed tie-breaks only on a score tie or unscored candidates (split read left for
    // supplementary-resolve). bwaHasPriority true leaves bwa's order untouched; mate is the other read's lifted record, or null.
    public static ApplyResult apply(
            final List<LiftedAlignment> alignments, final boolean concordant, final LiftedAlignment self,
            final int seed, final boolean bwaHasPriority, final LiftedRecord mate)
    {
        if(concordant || bwaHasPriority)
        {
            return keepBwaPrimary(alignments, self);
        }
        return pickByScore(alignments, self, seed, mate);
    }

    // Highest recomputed genome score wins ("score"). Top-score ties are settled in order by: mate proximity ("mate"),
    // then a spliced placement over a same-locus soft-clip ("junction"), then a read-name seed ("random"). Nothing is
    // dropped; losers ride in XA. Unscored candidates (a split read left for Step 3) keep bwa's primary.
    private static ApplyResult pickByScore(
            final List<LiftedAlignment> alignments, final LiftedAlignment self, final int seed, final LiftedRecord mate)
    {
        List<LiftedAlignment> candidates = new ArrayList<>();
        for(LiftedAlignment alignment : alignments)
        {
            if(!alignment.Dropped)
            {
                candidates.add(alignment);
            }
        }
        if(candidates.size() < 2)
        {
            return keepBwaPrimary(alignments, self);
        }

        int topScore = Integer.MIN_VALUE;
        for(LiftedAlignment alignment : candidates)
        {
            topScore = Math.max(topScore, alignment.GenomicScore);
        }
        if(topScore == Integer.MIN_VALUE)
        {
            return keepBwaPrimary(alignments, self); // unscored (split read left for Step 3): keep bwa's placement
        }

        // Collapse identical placements (same locus + CIGAR from different sources, e.g. a ref self and a tx alt that
        // lift to the same contiguous alignment) so the tie is over distinct placements, not weighted by source count.
        List<LiftedAlignment> top = new ArrayList<>();
        Set<String> topKeys = new HashSet<>();
        for(LiftedAlignment alignment : candidates)
        {
            if(alignment.GenomicScore == topScore && topKeys.add(placementKey(alignment)))
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
            List<LiftedAlignment> contenders = mateProximalSubset(top, mate);
            if(contenders.size() == 1)
            {
                winner = contenders.get(0);
                note = "mate";
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
        return new ApplyResult(indexOf(alignments, winner), winner, note);
    }

    private static ApplyResult keepBwaPrimary(final List<LiftedAlignment> alignments, final LiftedAlignment self)
    {
        return new ApplyResult(indexOf(alignments, self), self, "");
    }

    // identity, not equals: a candidate revised by withLiftedCigar is a distinct object at the same list position.
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

    // Tied candidates proximal to the mate; the full set unchanged when the mate is absent or does not discriminate.
    private static List<LiftedAlignment> mateProximalSubset(final List<LiftedAlignment> top, final LiftedRecord mate)
    {
        if(mate == null || !mate.hasPlacement())
        {
            return top;
        }
        List<LiftedAlignment> near = new ArrayList<>();
        for(LiftedAlignment alignment : top)
        {
            if(isMateProximal(alignment, mate))
            {
                near.add(alignment);
            }
        }
        return (near.isEmpty() || near.size() == top.size()) ? top : near;
    }

    private static boolean isMateProximal(final LiftedAlignment alignment, final LiftedRecord mate)
    {
        if(!mate.finalChromosome().equals(alignment.LiftedChromosome))
        {
            return false;
        }
        int mateEnd = mate.primaryAlignment().alignedEnd();
        int gap;
        if(alignment.LiftedPos > mateEnd)
        {
            gap = alignment.LiftedPos - mateEnd;
        }
        else if(alignment.LiftedPos < mate.finalPos())
        {
            gap = mate.finalPos() - alignment.LiftedPos;
        }
        else
        {
            gap = 0;
        }
        return gap <= MATE_PROXIMITY_MAX_DISTANCE;
    }

    // Tie-break within an equal-top-score set: a spliced placement (a real N junction) beats a clipped placement
    // (soft-clip, no N) at the same lifted locus. Both describe the same read at the same start; the junction is the
    // correct RNA interpretation - bwa soft-clipped rather than cross the intron - so it is not left to the coin.
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

    private static String placementKey(final LiftedAlignment alignment)
    {
        return alignment.locusKey() + ":" + alignment.LiftedCigar;
    }
}
