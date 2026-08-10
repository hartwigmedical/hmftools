package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MISMATCHES_AND_DELETIONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.XA_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.XS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.firstInPair;
import static com.hartwig.hmftools.tars.common.TarsConstants.SUPP_AS_DROP_THRESHOLD;
import static com.hartwig.hmftools.tars.liftback.SaTagRewriter.rewriteSaTag;

import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.tars.common.BwaScoring;
import com.hartwig.hmftools.tars.common.TarsCigarUtils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.SequenceUtil;

public final class BamRecordEmitter
{
    private static final String NUM_HITS_ATTRIBUTE = "NH";

    private final ContigTranslator mContigTranslator;
    private final boolean mSupplementaryResolverEnabled;
    private final RefGenomeInterface mRefGenome;
    private final ExcludedRegions mExcludedRegions;

    private int mPrimariesSeen;
    private int mPrimariesLiftFailed;

    public BamRecordEmitter(
            final ContigTranslator contigTranslator, final boolean supplementaryResolverEnabled,
            final RefGenomeInterface refGenome, final ExcludedRegions excludedRegions)
    {
        mContigTranslator = contigTranslator;
        mSupplementaryResolverEnabled = supplementaryResolverEnabled;
        mRefGenome = refGenome;
        mExcludedRegions = excludedRegions;
    }

    public int primariesSeen()
    {
        return mPrimariesSeen;
    }

    public int primariesLiftFailed()
    {
        return mPrimariesLiftFailed;
    }

    public void emit(
            final List<SAMRecord> records, final List<LiftedRecord> liftedRecords,
            final Set<Integer> absorbedSupplementaries, final boolean primaryUnmapped,
            final LiftedMatePair matePair, final Consumer<SAMRecord> consumer)
    {
        if(records.isEmpty())
        {
            return;
        }

        SAMRecord primary = records.get(0);
        LiftedRecord primaryResult = liftedRecords.get(0);
        boolean finalPrimaryUnmapped = primaryUnmapped || willBeUnmapped(primary, primaryResult);
        boolean[] willEmit = computeWillEmit(
                records, liftedRecords, absorbedSupplementaries, primary, finalPrimaryUnmapped);
        Set<String> droppedSuppSaKeys = SaTagRewriter.droppedSuppSaKeys(
                records, willEmit, primary, mContigTranslator);
        int numHits = LiftBackDiscriminator.countDistinctLoci(primaryResult);
        boolean primarySwapped = !finalPrimaryUnmapped && primaryResult.swapped();
        LiftedRecord[] resolved = liftedRecords.toArray(LiftedRecord[]::new);

        for(int i = 0; i < records.size(); ++i)
        {
            if(!willEmit[i])
            {
                continue;
            }

            SAMRecord record = records.get(i);
            String rebuiltSuppSa = primarySwapped && record != primary && record.getSupplementaryAlignmentFlag()
                    ? SaTagRewriter.buildSwappedSuppSa(records, resolved, willEmit, primaryResult, i)
                    : null;
            applyAndWrite(
                    record, liftedRecords.get(i), numHits, droppedSuppSaKeys, matePair, consumer,
                    rebuiltSuppSa, primaryUnmapped && record == primary);
        }
    }

    private boolean[] computeWillEmit(
            final List<SAMRecord> records, final List<LiftedRecord> liftedRecords,
            final Set<Integer> absorbedSupplementaries, final SAMRecord primary, final boolean primaryUnmapped)
    {
        Set<AlignmentKey> emittedSupplementaries = new HashSet<>();
        boolean[] willEmit = new boolean[records.size()];
        for(int i = 0; i < records.size(); ++i)
        {
            SAMRecord record = records.get(i);
            LiftedRecord result = liftedRecords.get(i);
            boolean drop = absorbedSupplementaries.contains(i);
            if(!drop && record.getSupplementaryAlignmentFlag() && primaryUnmapped)
            {
                drop = true;
            }
            if(!drop && record != primary && record.getSupplementaryAlignmentFlag()
                    && !emittedSupplementaries.add(dedupKey(result)))
            {
                drop = true;
            }
            if(!drop && mSupplementaryResolverEnabled && record.getSupplementaryAlignmentFlag()
                    && !record.getReadUnmappedFlag())
            {
                Integer alignmentScore = record.getIntegerAttribute(ALIGNMENT_SCORE_ATTRIBUTE);
                if(alignmentScore != null && alignmentScore < SUPP_AS_DROP_THRESHOLD)
                {
                    drop = true;
                }
            }
            if(!drop && record.getSupplementaryAlignmentFlag() && excludes(result))
            {
                drop = true;
            }
            willEmit[i] = !drop;
        }
        return willEmit;
    }

    private void applyAndWrite(
            final SAMRecord record, final LiftedRecord result, final int numHits,
            final Set<String> droppedSuppSaKeys, final LiftedMatePair matePair,
            final Consumer<SAMRecord> consumer, final String rebuiltSuppSa, final boolean unmapDecided)
    {
        boolean writeNumHits = result.hasPlacement() || !unmapDecided && !record.getReadUnmappedFlag();
        if(!record.isSecondaryOrSupplementary() && !record.getReadUnmappedFlag())
        {
            ++mPrimariesSeen;
            if(!result.hasPlacement() && !unmapDecided)
            {
                ++mPrimariesLiftFailed;
            }
        }

        applyResultToRecord(record, result, matePair, unmapDecided);
        if(record.getReadUnmappedFlag())
        {
            matePair.patchMateFields(record);
            if(!record.getSupplementaryAlignmentFlag())
            {
                consumer.accept(record);
            }
            return;
        }

        String rewrittenSa = rebuiltSuppSa != null
                ? rebuiltSuppSa
                : rewriteSaTag(record.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE), mContigTranslator, droppedSuppSaKeys);
        if(rewrittenSa == null && record.getSupplementaryAlignmentFlag())
        {
            return;
        }

        record.setAttribute(SUPPLEMENTARY_ATTRIBUTE, rewrittenSa);
        matePair.patchMateFields(record);
        if(writeNumHits)
        {
            record.setAttribute(NUM_HITS_ATTRIBUTE, numHits);
        }
        refreshNmDropMd(record);
        consumer.accept(record);
    }

    static void applyResultToRecord(
            final SAMRecord record, final LiftedRecord result, final LiftedMatePair matePair, final boolean unmapDecided)
    {
        if(!result.hasPlacement())
        {
            if(record.getReadUnmappedFlag())
            {
                return;
            }
            if(unmapDecided)
            {
                markPrimaryUnmapped(record, matePair);
                return;
            }
            if(record.isSecondaryOrSupplementary() && mirrorOwnPrimaryOntoFailedSupp(record, matePair))
            {
                return;
            }

            matePair.unmapRead(record);
            return;
        }

        Cigar liftedCigar = TarsCigarUtils.normalize(TextCigarCodec.decode(result.finalCigar()));
        if(result.negativeStrand() != record.getReadNegativeStrandFlag())
        {
            byte[] bases = record.getReadBases();
            if(bases != null && bases.length > 0)
            {
                byte[] flipped = bases.clone();
                Nucleotides.reverseComplementBasesInPlace(flipped, 0, flipped.length);
                record.setReadBases(flipped);
            }
            byte[] qualities = record.getBaseQualities();
            if(qualities != null && qualities.length > 0)
            {
                byte[] reversed = qualities.clone();
                for(int i = 0, j = reversed.length - 1; i < j; ++i, --j)
                {
                    byte value = reversed[i];
                    reversed[i] = reversed[j];
                    reversed[j] = value;
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
        record.setAttribute(XS_ATTRIBUTE, null);
        if(result.hasNCigar() && result.transcriptStrand() != 0)
        {
            record.setAttribute(XS_ATTRIBUTE, result.transcriptStrand() > 0 ? Character.valueOf('+') : Character.valueOf('-'));
        }
    }

    static void markPrimaryUnmapped(final SAMRecord record, final LiftedMatePair matePair)
    {
        matePair.unmapRead(record);
    }

    static boolean willBeUnmapped(final SAMRecord record, final LiftedRecord result)
    {
        return !result.hasPlacement() || exceedsMappingCap(record, result);
    }

    static boolean exceedsMappingCap(final SAMRecord record, final LiftedRecord result)
    {
        return record.getMappingQuality() == 0 && result.numXaAlts() == 0 && result.hasPlacement()
                && !result.primaryAlignment().FromTxContig;
    }

    public boolean excludes(final LiftedRecord result)
    {
        return mExcludedRegions != null && result.hasPlacement()
                && mExcludedRegions.excludes(
                        result.finalChromosome(), result.finalPos(), result.primaryAlignment().alignedEnd());
    }

    private static AlignmentKey dedupKey(final LiftedRecord result)
    {
        return result.hasPlacement()
                ? result.primaryAlignment().key()
                : new AlignmentKey("*", 0, "*", true);
    }

    private static boolean mirrorOwnPrimaryOntoFailedSupp(final SAMRecord record, final LiftedMatePair matePair)
    {
        if(!record.getReadPairedFlag())
        {
            return false;
        }
        LiftedRecord primary = matePair.ownPrimary(firstInPair(record));
        if(primary == null || !primary.hasPlacement())
        {
            return false;
        }

        record.setReferenceName(primary.finalChromosome());
        record.setAlignmentStart(primary.finalPos());
        record.setCigarString(primary.finalCigar());
        record.setMappingQuality(0);
        record.setAttribute(XA_ATTRIBUTE, null);
        return true;
    }

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

    private int computeEditDistance(final SAMRecord record)
    {
        byte[] readBases = record.getReadBases();
        if(readBases == null || readBases.length == 0)
        {
            return -1;
        }

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
                    if(readIndex + length > readBases.length)
                    {
                        return -1;
                    }
                    byte[] refBases = BwaScoring.refWindow(
                            mRefGenome, record.getReferenceName(), refPos, refPos + length - 1);
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
