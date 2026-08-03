package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MISMATCHES_AND_DELETIONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.XA_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.XS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.firstInPair;
import static com.hartwig.hmftools.tars.common.TarsConstants.CONFIDENT_MAPQ;
import static com.hartwig.hmftools.tars.common.TarsConstants.PRIMARY_AS_UNMAP_THRESHOLD;
import static com.hartwig.hmftools.tars.common.TarsConstants.SUPP_AS_DROP_THRESHOLD;
import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;
import static com.hartwig.hmftools.tars.liftback.SaTagRewriter.rewriteSaTag;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.common.BwaScoring;
import com.hartwig.hmftools.tars.liftback.overhang.OverhangGate;
import com.hartwig.hmftools.tars.liftback.supplementary.SupplementaryResolver;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.SequenceUtil;

//TODO: revisit this again

// Per-name-group transform engine: lifts one contiguous read-name group to genomic coordinates, runs the optional
// supplementary-merge and overhang-gate passes, patches mate fields, and hands each kept record to the consumer.
// A name-sorted group holds both mates and all their supplementaries, so the per-group LiftedMatePair is self-sufficient.
public class LiftBackGroupProcessor
{
    // hmf-common has no NH constant and this file holds the only two uses in tars
    private static final String NUM_HITS_ATTRIBUTE = "NH";

    private final LiftBackDiscriminator mDiscriminator;

    // optional supplementary resolver: merges primary + supp across annotated junctions. Null when disabled.
    private final SupplementaryResolver mSupplementaryResolver;

    // wraps the overhang gate over a record's lifted alignments (collapse candidates, finalize chosen primary + alts).
    private final OverhangReconciler mOverhangReconciler;

    // genomic reference for the post-lift NM recompute; null when no reference is loaded, in which case NM is cleared.
    private final RefGenomeInterface mRefGenome;

    // genomic excluded regions (rRNA / contamination), checked post-lift. Null when no regions file configured.
    private final ExcludedRegions mExcludedRegions;

    private int mRecordsSeen;

    // primaries that arrived mapped and end up unmapped after lift-back; a high rate means a sidecar/FASTA mismatch.
    private int mPrimariesUnmapped;

    public LiftBackGroupProcessor(
            final LiftBackDiscriminator discriminator, final SupplementaryResolver supplementaryResolver,
            final OverhangGate overhangGate, final RefGenomeInterface refGenome, final ExcludedRegions excludedRegions)
    {
        mDiscriminator = discriminator;
        mSupplementaryResolver = supplementaryResolver;
        mOverhangReconciler = new OverhangReconciler(overhangGate);
        mRefGenome = refGenome;
        mExcludedRegions = excludedRegions;
    }

    public int recordsSeen() { return mRecordsSeen; }

    public int primariesUnmapped() { return mPrimariesUnmapped; }

    // Processes all records sharing one read name as a group so primary + supps resolve together. This lets
    // supplementary resolve see all split-read components and lets read /2 hint read /1's resolved introns.
    public void processNameGroup(final List<SAMRecord> group, final Consumer<SAMRecord> consumer)
    {
        mRecordsSeen += group.size();

        LiftedMatePair matePair = new LiftedMatePair();
        List<SAMRecord> firstOfPair = new ArrayList<>();
        List<SAMRecord> secondOfPair = new ArrayList<>();
        for(SAMRecord record : group)
        {
            if(firstInPair(record))
            {
                firstOfPair.add(record);
            }
            else
            {
                secondOfPair.add(record);
            }
        }

        // Decide /1 first, then /2 with /1's resolved junctions as intron hints and /1's lifted locus for
        // mate-proximity. Hinting is one-way: /2 sees /1; /1 is not re-decided against /2.
        // /1 is decided before /2 has been lifted, so it has no mate to be placed against.
        ReadDecision firstOfPairDecision = decideRead(firstOfPair, Collections.emptyList(), null);
        recordPrimaryInPair(firstOfPair, firstOfPairDecision, matePair);
        ReadDecision secondOfPairDecision = decideRead(secondOfPair, firstOfPairDecision.IntroducedIntrons, matePair.mateOf(false));
        recordPrimaryInPair(secondOfPair, secondOfPairDecision, matePair);

        writeRead(firstOfPair, firstOfPairDecision, matePair, consumer);
        writeRead(secondOfPair, secondOfPairDecision, matePair, consumer);
    }

    // Records the post-supplementary-resolve primary so the mate's MC tag uses the merged/extended cigar rather
    // than the pre-resolve one the seed lift produced.
    private static void recordPrimaryInPair(
            final List<SAMRecord> readRecords, final ReadDecision decision, final LiftedMatePair pair)
    {
        LiftedRecord primaryResult = decision.primaryResult();
        if(readRecords.isEmpty() || primaryResult == null)
        {
            return;
        }
        pair.recordPrimary(firstInPair(readRecords.get(0)), primaryResult);
    }

    // Per-read decision: resolves every record, runs the supplementary resolver, returns the lifted results +
    // per-record drop flags + introns introduced by supplementary resolve (for mate hinting). Does not emit.
    private static final class ReadDecision
    {
        LiftedRecord[] LiftedRecords;
        boolean[] DroppedBySupplementary;
        List<ChrBaseRegion> IntroducedIntrons;

        // set when decide flipped the primary to unmapped, which the emit step treats differently from a failed lift
        boolean PrimaryUnmapped;

        ReadDecision(final LiftedRecord[] liftedRecords, final boolean[] droppedBySupplementary,
                final List<ChrBaseRegion> introducedIntrons)
        {
            LiftedRecords = liftedRecords;
            DroppedBySupplementary = droppedBySupplementary;
            IntroducedIntrons = introducedIntrons != null ? introducedIntrons : Collections.emptyList();
        }

        LiftedRecord primaryResult()
        {
            return LiftedRecords.length > 0 ? LiftedRecords[0] : null;
        }
    }

    private ReadDecision decideRead(
            final List<SAMRecord> records, final List<ChrBaseRegion> mateHintIntrons, final LiftedRecord mate)
    {
        if(records.isEmpty())
        {
            return new ReadDecision(new LiftedRecord[0], null, Collections.emptyList());
        }

        // primary is the first record (bwa emits a read's primary before its supplementaries); crash otherwise.
        SAMRecord primary = records.get(0);
        if(primary.getSupplementaryAlignmentFlag())
        {
            throw new IllegalStateException(String.format(
                    "read %s: first record in group is a supplementary, expected the primary first",
                    primary.getReadName()));
        }

        LiftedRecord primaryResult = mDiscriminator.resolve(primary, mOverhangReconciler, mate);

        LiftedRecord[] resolved = new LiftedRecord[records.size()];
        for(int i = 0; i < records.size(); ++i)
        {
            SAMRecord record = records.get(i);
            resolved[i] = record == primary ? primaryResult : mDiscriminator.resolve(record, mOverhangReconciler);
        }

        List<ChrBaseRegion> introduced = new ArrayList<>();
        boolean[] droppedBySupplementary = applySupplementaryResolve(records, resolved, primary, mateHintIntrons, introduced);

        // Finalize the chosen primary: collapse, then extend the terminal softclip (tx-match only). Does real work beyond the
        // per-candidate collapse: (1) supplementary resolve can fabricate a fresh M N M with its own terminal micro-junction;
        // (2) the softclip extension is deferred here (post-resolve) so resolve saw the original clip. Clean primary unchanged.
        mOverhangReconciler.reconcileChosenPrimary(resolved, primary, 0, droppedBySupplementary);

        ReadDecision decision = new ReadDecision(resolved, droppedBySupplementary, introduced);

        // post-lift exclusion: if the final primary placement lands in an excluded (rRNA / contamination) region,
        // unmap it REDUX-style - flip the result to unplaced so the mate is coordinated through the pair (willBeUnmapped)
        // and the record is unmapped at apply time. Lifted genomic coords are the only space tx-contig reads can be
        // tested against the genomic region list, so this is post-lift, not the old pre-lift fragment pre-filter.
        if(liftsIntoExcludedRegion(resolved[0]))
        {
            resolved[0] = LiftedRecord.unmapped("excluded_region_unmapped");
            decision.PrimaryUnmapped = true;
        }
        else if(!primary.getReadUnmappedFlag())
        {
            // Unmap the chosen primary during decide (not emit) so the mate pair and supp-orphan gate observe it. Two
            // triggers: over bwa's XA cap (MAPQ 0, no XA = too many loci); or a residual short anchor below bwa's -T 30
            // floor liftback couldn't improve, only when supp resolve left the primary's stale bwa AS valid.
            boolean primaryPostProcessed = resolved[0] != primaryResult;
            if(exceedsMappingCap(primary, resolved[0]))
            {
                resolved[0] = LiftedRecord.unmapped("over_cap_unmapped");
                decision.PrimaryUnmapped = true;
                TARS_LOGGER.trace("over-cap unmap {}: inputMapQuality=0, no XA (past bwa XA cap)", primary.getReadName());
            }
            else if(mSupplementaryResolver != null && !primaryPostProcessed)
            {
                Integer alignmentScore = primary.getIntegerAttribute(ALIGNMENT_SCORE_ATTRIBUTE);
                if(alignmentScore != null && alignmentScore < PRIMARY_AS_UNMAP_THRESHOLD)
                {
                    resolved[0] = LiftedRecord.unmapped("low_as_unmapped");
                    decision.PrimaryUnmapped = true;
                    TARS_LOGGER.trace(
                            "AS-floor unmap {}: AS={} < {}", primary.getReadName(), alignmentScore, PRIMARY_AS_UNMAP_THRESHOLD);
                }
            }
        }

        return decision;
    }

    private void writeRead(
            final List<SAMRecord> records, final ReadDecision decision, final LiftedMatePair matePair,
            final Consumer<SAMRecord> consumer)
    {
        if(records.isEmpty())
        {
            return;
        }

        LiftedRecord[] resolved = decision.LiftedRecords;
        boolean[] droppedBySupplementary = decision.DroppedBySupplementary;
        LiftedRecord primaryResult = decision.primaryResult();
        SAMRecord primary = records.get(0);

        // A supplementary whose primary ends up unmapped (e.g. the primary lifted into an excluded region) is
        // orphaned: there is no emitted primary for its SA to reference, so REDUX FragmentCoords would read a
        // dangling primary locus off it. Drop the whole group's supps in that case.
        boolean primaryUnmapped = decision.PrimaryUnmapped || willBeUnmapped(primary, primaryResult);

        boolean[] willEmit = computeWillEmit(records, resolved, droppedBySupplementary, primary, primaryUnmapped);

        // Remove dropped supps' SA entries from the primary so it never references a missing supp.
        Set<String> droppedSuppSaKeys = SaTagRewriter.droppedSuppSaKeys(
                records, willEmit, primary, mDiscriminator.contigTranslator());

        // NH = distinct genomic loci, recomputed here because reconcileChosenPrimary's alt-softclip extension can collapse
        // an alt onto the primary after numLoci was baked. Recounting keeps NH consistent with the XA the record carries;
        // counting emitted records would inflate it by tx-contig multiplicity.
        int nh = LiftBackDiscriminator.countDistinctLoci(primaryResult);

        // When the discriminator moved the primary off bwa's placement (a ref->tx swap or a random-locus pick), each
        // surviving supplementary's own bwa SA still references bwa's PRE-move primary locus, which was never emitted.
        // Rebuild those supp SAs to reference the emitted primary so REDUX FragmentCoords reads the right primary coords.
        boolean primarySwapped = !primaryUnmapped && primaryResult.swapped();

        for(int i = 0; i < records.size(); ++i)
        {
            if(!willEmit[i])
            {
                continue;
            }

            SAMRecord record = records.get(i);
            String rebuiltSuppSa = (primarySwapped && record != primary && record.getSupplementaryAlignmentFlag())
                    ? SaTagRewriter.buildSwappedSuppSa(records, resolved, willEmit, primaryResult, i) : null;
            applyAndWriteRecord(
                    record, resolved[i], nh, droppedSuppSaKeys, matePair, consumer, rebuiltSuppSa,
                    decision.PrimaryUnmapped && record == primary);
        }
    }

    // Decides which records in the group are emitted vs dropped. Drops: absorbed-by-supp-merge, orphaned supp (primary
    // unmapped), duplicate supp (same lifted locus across tx contigs), low-AS -T-19 supp, and excluded-region supp.
    private boolean[] computeWillEmit(
            final List<SAMRecord> records, final LiftedRecord[] resolved, final boolean[] droppedBySupplementary,
            final SAMRecord primary, final boolean primaryUnmapped)
    {
        Set<String> emittedSuppKeys = new HashSet<>();
        boolean[] willEmit = new boolean[records.size()];
        for(int i = 0; i < records.size(); ++i)
        {
            SAMRecord record = records.get(i);
            LiftedRecord result = resolved[i];
            boolean drop = droppedBySupplementary != null && droppedBySupplementary[i];
            if(!drop && record.getSupplementaryAlignmentFlag() && primaryUnmapped)
            {
                drop = true;
            }
            if(!drop && record != primary && record.getSupplementaryAlignmentFlag())
            {
                String key = dedupKey(result);
                if(!emittedSuppKeys.add(key))
                {
                    drop = true;
                }
            }
            // Drop supps the supplementary-resolve pass left behind that exist only because bwa-mem2 was run with
            // -T 19 below its default of 30 (see SUPP_AS_DROP_THRESHOLD). Only applies when supplementary resolve
            // ran, because a configuration without it might want to retain these supps for other reasons.
            if(!drop && mSupplementaryResolver != null && record.getSupplementaryAlignmentFlag()
                    && !record.getReadUnmappedFlag())
            {
                Integer alignmentScore = record.getIntegerAttribute(ALIGNMENT_SCORE_ATTRIBUTE);
                if(alignmentScore != null && alignmentScore < SUPP_AS_DROP_THRESHOLD)
                {
                    drop = true;
                }
            }
            // a supp lifting into an excluded region is dropped (a supp can't be unmapped); its SA entry is
            // removed from the primary below so the primary doesn't reference a supp that isn't emitted.
            if(!drop && record.getSupplementaryAlignmentFlag() && liftsIntoExcludedRegion(result))
            {
                drop = true;
            }
            willEmit[i] = !drop;
        }
        return willEmit;
    }

    // Builds a supplementary-merge candidate from the post-lift primary + supps and applies the merge if accepted. Returns a
    // parallel array marking records absorbed by the merge (drop from emission), or null when disabled / no-op.
    // introducedIntronsOut is appended-to so the partner mate can use the introns as hints in the second pass.
    private boolean[] applySupplementaryResolve(
            final List<SAMRecord> records, final LiftedRecord[] resolved, final SAMRecord primary,
            final List<ChrBaseRegion> mateHintIntrons, final List<ChrBaseRegion> introducedIntronsOut)
    {
        if(mSupplementaryResolver == null || primary.getReadUnmappedFlag())
        {
            return null;
        }

        int primaryIdx = -1;
        List<Integer> suppIndices = new ArrayList<>();
        for(int i = 0; i < records.size(); ++i)
        {
            SAMRecord record = records.get(i);
            if(record == primary)
            {
                primaryIdx = i;
            }
            else if(record.getSupplementaryAlignmentFlag() && !record.getReadUnmappedFlag())
            {
                suppIndices.add(i);
            }
        }
        if(primaryIdx < 0)
        {
            return null;
        }

        LiftedRecord primaryRes = resolved[primaryIdx];
        if(!primaryRes.hasPlacement())
        {
            return null;
        }

        List<SupplementaryResolver.Supplementary> suppDtos = new ArrayList<>(suppIndices.size());
        for(int i = 0; i < suppIndices.size(); ++i)
        {
            int suppRecordIndex = suppIndices.get(i);
            SAMRecord supp = records.get(suppRecordIndex);
            LiftedRecord suppRes = resolved[suppRecordIndex];
            if(!suppRes.hasPlacement())
            {
                continue;
            }
            suppDtos.add(new SupplementaryResolver.Supplementary(
                    i, suppRes.finalChromosome(), !supp.getReadNegativeStrandFlag(),
                    suppRes.finalPos(), suppRes.finalCigar(), supp.getMappingQuality()));
        }

        SupplementaryResolver.Candidate candidate = new SupplementaryResolver.Candidate(
                primaryRes.finalChromosome(), !primary.getReadNegativeStrandFlag(), primary.getReadLength(),
                primaryRes.finalPos(), primaryRes.finalCigar(), suppDtos,
                primary.getReadBases(), mateHintIntrons);

        SupplementaryResolver.Result result = mSupplementaryResolver.resolve(candidate);
        if(!result.merged())
        {
            return null;
        }

        // Rewrite the primary's result with the merged (spliced) cigar + start (start may shift on a left-extend).
        // Merged MAPQ = max of primary and every merged supp, bumped to 60 only when the pair maps to a single locus
        // (primary lifts to one locus, no competing XA alt); a pair at >1 locus keeps its bwa MAPQ.
        int mergedMapQuality = primaryRes.updatedMapQuality();
        for(Integer dtoIdx : result.droppedSupplementaryIndices())
        {
            mergedMapQuality = Math.max(mergedMapQuality, records.get(suppIndices.get(dtoIdx)).getMappingQuality());
        }
        if(primaryRes.numLoci() == 1)
        {
            mergedMapQuality = CONFIDENT_MAPQ;
        }
        resolved[primaryIdx] = primaryRes.withRevisedPrimary(
                result.mergedStart(), result.mergedCigar(), mergedMapQuality, "supplementary-resolved");
        TARS_LOGGER.trace(
                "supplementary-resolved {}: {}:{} {} -> {} (depth {})",
                primary.getReadName(), primaryRes.finalChromosome(), primaryRes.finalPos(),
                primaryRes.finalCigar(), result.mergedCigar(), result.chainDepth());

        if(introducedIntronsOut != null)
        {
            introducedIntronsOut.addAll(result.introducedIntrons());
        }

        boolean[] dropped = new boolean[records.size()];
        for(Integer dtoIdx : result.droppedSupplementaryIndices())
        {
            dropped[suppIndices.get(dtoIdx)] = true;
        }
        return dropped;
    }

    // Dedup key for supplementaries: lifted (chrom, pos, cigar) + lifted strand. Opposite-strand placements at the
    // same coords/cigar are kept as distinct records. Every unplaced supplementary keys the same, so only the first
    // of them survives the dedup.
    private static String dedupKey(final LiftedRecord result)
    {
        if(!result.hasPlacement())
        {
            return "*:0:*:+";
        }

        return result.finalChromosome() + ":" + result.finalPos() + ":" + result.finalCigar()
                + ":" + (result.negativeStrand() ? '-' : '+');
    }

    // true if the result's lifted genomic span overlaps an excluded region. False for unplaced results.
    private boolean liftsIntoExcludedRegion(final LiftedRecord result)
    {
        if(mExcludedRegions == null || !result.hasPlacement())
        {
            return false;
        }

        return mExcludedRegions.excludes(result.finalChromosome(), result.finalPos(), result.primaryAlignment().alignedEnd());
    }

    private void applyAndWriteRecord(
            final SAMRecord record, final LiftedRecord result, final int nh, final Set<String> droppedSuppSaKeys,
            final LiftedMatePair matePair, final Consumer<SAMRecord> consumer, final String rebuiltSuppSa,
            final boolean unmapDecided)
    {
        // Captured before the record is mutated: a record whose lift failed is still emitted and still carries NH,
        // only a record that arrived unmapped or was unmapped at decide time goes without.
        boolean writeNumHits = result.hasPlacement() || !unmapDecided && !record.getReadUnmappedFlag();

        // Also pre-mutation, and the single place both unmap routes pass through: a failed lift and a decide-time flip.
        if(!result.hasPlacement() && !record.isSecondaryOrSupplementary() && !record.getReadUnmappedFlag())
        {
            ++mPrimariesUnmapped;
        }

        applyResultToRecord(record, result, matePair, unmapDecided);

        if(record.getReadUnmappedFlag())
        {
            record.setReadNegativeStrandFlag(false);
        }

        // rebuiltSuppSa is set only for a supplementary whose primary the discriminator moved: use the primary-referencing
        // SA rebuilt from the emitted group. Otherwise translate this record's own bwa SA to genomic coords as usual.
        String rewrittenSa = rebuiltSuppSa != null
                ? rebuiltSuppSa
                : rewriteSaTag(
                        record.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE), mDiscriminator.contigTranslator(), droppedSuppSaKeys);

        // A record that stays supplementary must carry an SA tag: REDUX dedup (FragmentCoords.fromRead)
        // reads the primary's coords from it. If every SA entry failed to lift the supplementary is orphaned
        // in genomic space and cannot be represented, so drop it rather than emit a supp with a null SA.
        if(rewrittenSa == null && record.getSupplementaryAlignmentFlag())
        {
            return;
        }

        record.setAttribute(SUPPLEMENTARY_ATTRIBUTE, rewrittenSa);
        matePair.patchMateFields(record);

        if(writeNumHits)
        {
            record.setAttribute(NUM_HITS_ATTRIBUTE, nh);
        }

        refreshNmDropMd(record);

        consumer.accept(record);
    }

    // Supplementaries whose own lift failed are mirrored onto their primary's lifted coords (keeping the 0x800 flag)
    // rather than marked unmapped, since htsjdk's validator rejects 0x4+0x800.
    private static void applyResultToRecord(
            final SAMRecord record, final LiftedRecord result, final LiftedMatePair matePair, final boolean unmapDecided)
    {
        if(!result.hasPlacement())
        {
            // arrived unmapped: nothing was lifted and nothing needs clearing
            if(record.getReadUnmappedFlag())
            {
                return;
            }

            // flipped to unmapped at decide time (excluded region / over-cap / low-AS): a fuller clear than a failed
            // lift, since the record still carries coords, NH and mate state that the flip invalidates.
            if(unmapDecided)
            {
                markPrimaryUnmapped(record);
                return;
            }

            if(record.isSecondaryOrSupplementary() && mirrorOwnPrimaryOntoFailedSupp(record, matePair))
            {
                return;
            }

            record.setReadUnmappedFlag(true);
            record.setReadNegativeStrandFlag(false);
            record.setReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
            record.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
            record.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
            record.setMappingQuality(0);
            record.setAttribute(XA_ATTRIBUTE, null);
            record.setAttribute(SUPPLEMENTARY_ATTRIBUTE, null);
            return;
        }

        Cigar liftedCigar = new Cigar(ContigTranslator.mergeAdjacentSameOp(
                ContigTranslator.dropZeroLength(TextCigarCodec.decode(result.finalCigar()).getCigarElements())));

        // A discriminator swap can flip the primary onto an opposite-strand placement (bwa's contiguous
        // alignment vs the winning spliced tx placement). Reverse-complement bases and reverse quals so SEQ stays
        // on the genomic forward strand matching the lifted pos/cigar; else refreshNmDropMd's NM inflates to ~read length.
        if(result.negativeStrand() != record.getReadNegativeStrandFlag())
        {
            byte[] bases = record.getReadBases();
            if(bases != null && bases.length > 0)
            {
                byte[] flipped = bases.clone();
                Nucleotides.reverseComplementBasesInPlace(flipped, 0, flipped.length);
                record.setReadBases(flipped);
            }
            byte[] quals = record.getBaseQualities();
            if(quals != null && quals.length > 0)
            {
                byte[] reversed = quals.clone();
                for(int i = 0, j = reversed.length - 1; i < j; ++i, --j)
                {
                    byte swapValue = reversed[i];
                    reversed[i] = reversed[j];
                    reversed[j] = swapValue;
                }
                record.setBaseQualities(reversed);
            }
        }

        record.setReferenceName(result.finalChromosome());
        record.setAlignmentStart(result.finalPos());
        record.setCigar(liftedCigar);
        record.setReadNegativeStrandFlag(result.negativeStrand());
        record.setMappingQuality(result.updatedMapQuality());
        record.setAttribute(XA_ATTRIBUTE, result.xaTag());
        // XS:A:+/- carries transcript strand for Isofox's stranded junction interpretation. Clear first because
        // bwa-mem2 reuses the XS tag for its XS:i: sub-optimal score. Set only when the lifted cigar has an N AND
        // the transcript strand is known (tx-contig primary); ref-only N-cigars have no strand threaded, ship without XS.
        record.setAttribute(XS_ATTRIBUTE, null);
        if(result.hasNCigar() && result.transcriptStrand() != 0)
        {
            record.setAttribute(XS_ATTRIBUTE, result.transcriptStrand() > 0 ? Character.valueOf('+') : Character.valueOf('-'));
        }
    }

    // Paired-only: single-end bwa runs without -Y, so its supplementaries are hard-clipped and grafting the primary's
    // full-read cigar onto one would make computeEditDistance over-read the read bases. Those drop as orphans instead.
    private static boolean mirrorOwnPrimaryOntoFailedSupp(final SAMRecord record, final LiftedMatePair matePair)
    {
        if(!record.getReadPairedFlag())
        {
            return false;
        }

        LiftedRecord ownPrimary = matePair.ownPrimary(firstInPair(record));
        if(ownPrimary == null || !ownPrimary.hasPlacement())
        {
            return false;
        }

        record.setReferenceName(ownPrimary.finalChromosome());
        record.setAlignmentStart(ownPrimary.finalPos());
        record.setCigarString(ownPrimary.finalCigar());
        record.setMappingQuality(0);
        record.setAttribute(XA_ATTRIBUTE, null);
        // SA stays: rewriteSaTag translates tx-contig entries downstream and FragmentCoords requires it on supps
        return true;
    }

    // Mark a previously-mapped primary as unmapped, clearing every per-record tag the SAM spec invariant
    // would otherwise leave inconsistent. SA/XA/NH/MC all reference now-stale coords; ProperPair and
    // TLEN are meaningless once one end is unmapped.
    private static void markPrimaryUnmapped(final SAMRecord record)
    {
        record.setReadUnmappedFlag(true);
        record.setReadNegativeStrandFlag(false);
        record.setReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
        record.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
        record.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
        record.setMappingQuality(0);
        record.setAttribute(SUPPLEMENTARY_ATTRIBUTE, null);
        record.setAttribute(XA_ATTRIBUTE, null);
        record.setAttribute(NUM_HITS_ATTRIBUTE, null);
        record.setAttribute(MATE_CIGAR_ATTRIBUTE, null);
        if(record.getReadPairedFlag())
        {
            record.setProperPairFlag(false);
            record.setInferredInsertSize(0);
        }
    }

    // Single source of truth for whether a primary's final state ends up unmapped, so the mate patch and the
    // supp-orphan gate agree with what the record itself is about to become.
    private static boolean willBeUnmapped(final SAMRecord record, final LiftedRecord result)
    {
        return !result.hasPlacement() || exceedsMappingCap(record, result);
    }

    // bwa caps the XA list (-h 75); an over-cap read is emitted MAPQ 0 with NO XA (suppressed, not truncated), too many
    // places to trust, so unmap it. Input MAPQ 0 + no lifted alt distinguishes it from a unique read and a few-way multimapper.
    // A tx-contig primary is carved out: its 75+ hits are one gene's shared-exon contigs lifting to one locus, not distinct loci.
    private static boolean exceedsMappingCap(final SAMRecord record, final LiftedRecord result)
    {
        return record.getMappingQuality() == 0 && result.numXaAlts() == 0 && result.hasPlacement()
                && !result.primaryAlignment().FromTxContig;
    }

    // tx-contig MD/NM are stale once the read is lifted and recut. MD is dropped not rebuilt (rebuilding across a spliced
    // read needs reference spanning the whole intron); NM is recomputed cheaply from the aligned M blocks (N/S/H contribute
    // nothing) against the genomic reference. With no reference loaded or an unmapped record, NM is cleared.
    private void refreshNmDropMd(final SAMRecord record)
    {
        record.setAttribute(MISMATCHES_AND_DELETIONS_ATTRIBUTE, null);

        if(record.getReadUnmappedFlag() || mRefGenome == null)
        {
            record.setAttribute(NUM_MUTATONS_ATTRIBUTE, null);
            return;
        }

        int editDistance = computeEditDistance(record);
        record.setAttribute(NUM_MUTATONS_ATTRIBUTE, editDistance >= 0 ? Integer.valueOf(editDistance) : null);
    }

    // NM per the SAM spec: mismatches in aligned blocks plus inserted and deleted bases. N/S/H/P contribute nothing, so
    // each M block's reference is fetched on its own and the intron is never read. Returns -1 if any block's reference
    // is unavailable, so the caller leaves NM unset rather than write a wrong value.
    private int computeEditDistance(final SAMRecord record)
    {
        byte[] readBases = record.getReadBases();
        if(readBases == null || readBases.length == 0)
        {
            return -1;
        }

        String chromosome = record.getReferenceName();
        int refPos = record.getAlignmentStart();
        int readIndex = 0;
        int editDistance = 0;

        for(CigarElement element : record.getCigar().getCigarElements())
        {
            int length = element.getLength();
            switch(element.getOperator())
            {
                case M:
                case EQ:
                case X:
                    // Guard against a cigar whose read span exceeds the record's SEQ length (e.g. a failed-supp
                    // mirror that took its primary's full-length cigar onto a hard-clipped, shorter SEQ). Leave NM
                    // unset rather than index past the read bases.
                    if(readIndex + length > readBases.length)
                    {
                        return -1;
                    }
                    byte[] refBases = BwaScoring.refWindow(mRefGenome, chromosome, refPos, refPos + length - 1);
                    if(refBases == null)
                    {
                        return -1;
                    }
                    for(int i = 0; i < length; ++i)
                    {
                        if(!SequenceUtil.basesEqual(readBases[readIndex + i], refBases[i]))
                        {
                            ++editDistance;
                        }
                    }
                    refPos += length;
                    readIndex += length;
                    break;

                case I:
                    editDistance += length;
                    readIndex += length;
                    break;

                case D:
                    editDistance += length;
                    refPos += length;
                    break;

                case N:
                    refPos += length;
                    break;

                case S:
                    readIndex += length;
                    break;

                case H:
                case P:
                    break;
            }
        }

        return editDistance;
    }
}
