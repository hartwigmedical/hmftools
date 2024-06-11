package com.hartwig.hmftools.sage.sync;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.sync.FragmentSyncType.CIGAR_MISMATCH;
import static com.hartwig.hmftools.sage.sync.FragmentSyncType.NO_OVERLAP;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public class FragmentSync
{
    private final FragmentSyncReadHandler mReadHandler;
    private final Map<String,SAMRecord> mCachedReads;
    private final SAMFileHeader mSAMHeader;

    private final int[] mSyncCounts;

    public FragmentSync(final FragmentSyncReadHandler readHandler, final RefGenomeInterface refGenomeSource)
    {
        mReadHandler = readHandler;
        mCachedReads = Maps.newHashMap();
        mSyncCounts = new int[FragmentSyncType.values().length];
        mSAMHeader = buildRefSamHeader(refGenomeSource);
    }

    public void emptyCachedReads()
    {
        if(!mCachedReads.isEmpty())
        {
            mCachedReads.values().forEach(x -> mReadHandler.processReadRecord(x, false, null));
            mCachedReads.clear();
        }
    }

    public final int[] getSynCounts() { return mSyncCounts; }

    public boolean handleOverlappingReads(final SAMRecord record)
    {
        if(!record.getReadPairedFlag() || record.getMateUnmappedFlag())
            return false;

        final SAMRecord otherRecord = mCachedReads.get(record.getReadName());

        if(otherRecord != null)
        {
            try
            {
                FragmentSyncOutcome syncOutcome = CombinedSyncData.formFragmentRead(otherRecord, record, mSAMHeader);
                mCachedReads.remove(record.getReadName());

                SAMRecord fragmentRecord = syncOutcome.CombinedRecord;
                ++mSyncCounts[syncOutcome.SyncType.ordinal()];

                /*
                SG_LOGGER.trace("fragment sync: first({} {}:{}-{} {}) second({} {}:{}-{} {}) outcome({}) {}",
                        otherRecord.getReadName(), otherRecord.getContig(), otherRecord.getAlignmentStart(), otherRecord.getAlignmentEnd(),
                        otherRecord.getCigarString(), record.getReadName(), record.getContig(), record.getAlignmentStart(),
                        record.getAlignmentEnd(), record.getCigarString(), syncOutcome.SyncType,
                        syncOutcome.CombinedRecord != null ? format("newCigar(%s)", syncOutcome.CombinedRecord.getCigarString()) : "");
                */

                if(fragmentRecord != null)
                {
                    FragmentData fragmentData = new FragmentData(otherRecord, record);

                    mReadHandler.processReadRecord(fragmentRecord, false, fragmentData);
                }
                else if(syncOutcome.SyncType == CIGAR_MISMATCH)
                {
                    // favour the read with the longest INDEL where they disagree
                    int firstIndelLen = otherRecord.getCigar().getCigarElements().stream()
                            .filter(x -> x.getOperator().isIndel()).mapToInt(x -> x.getLength()).sum();

                    int secondIndelLen = record.getCigar().getCigarElements().stream()
                            .filter(x -> x.getOperator().isIndel()).mapToInt(x -> x.getLength()).sum();

                    if(secondIndelLen > firstIndelLen)
                        mReadHandler.processReadRecord(record, false);
                    else
                        mReadHandler.processReadRecord(otherRecord, false);
                }
                else if(syncOutcome.SyncType.processSeparately())
                {
                    mReadHandler.processReadRecord(otherRecord, false);
                    mReadHandler.processReadRecord(record, false);
                }
                else
                {
                    // only the first record
                    mReadHandler.processReadRecord(otherRecord, false);
                }
            }
            catch(Exception e)
            {
                ++mSyncCounts[FragmentSyncType.EXCEPTION.ordinal()];

                SG_LOGGER.error("failed to sync fragments: {}", e.toString());
                e.printStackTrace();

                SG_LOGGER.info("firstRead({} {}:{}-{} {})",
                        otherRecord.getReadName(), otherRecord.getContig(), otherRecord.getAlignmentStart(), otherRecord.getAlignmentEnd(),
                        otherRecord.getCigarString());

                SG_LOGGER.info("secondRead({} {}:{}-{} {})",
                        record.getReadName(), record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd(),
                        record.getCigarString());
            }

            return true;
        }

        // no cache for reads where the mate doesn't overlap
        if(!record.getContig().equals(record.getMateReferenceName()))
            return false;

        // if the mate is earlier, then it should have been processed and so no point in not handling this read now
        if(record.getMateAlignmentStart() < record.getAlignmentStart())
            return false;

        if(!positionsOverlap(
                record.getAlignmentStart(), record.getAlignmentEnd(),
                record.getMateAlignmentStart(), record.getMateAlignmentStart() + record.getReadLength()))
        {
            ++mSyncCounts[NO_OVERLAP.ordinal()];
            return false;
        }

        // cache until the paired read arrives
        mCachedReads.put(record.getReadName(), record);
        return true;
    }

    private SAMFileHeader buildRefSamHeader(final RefGenomeInterface refGenomeInterface)
    {
        if(!(refGenomeInterface instanceof RefGenomeSource))
            return null;

        RefGenomeSource refGenomeSource = (RefGenomeSource)refGenomeInterface;
        List<SAMSequenceRecord> sequenceRecords = Lists.newArrayList();

        SAMSequenceDictionary refGenomeDict = refGenomeSource.refGenomeFile().getSequenceDictionary();
        SAMSequenceDictionary dictionary = new SAMSequenceDictionary();

        for(SAMSequenceRecord sequence : refGenomeDict.getSequences())
        {
            if(HumanChromosome.contains(sequence.getSequenceName()) || MitochondrialChromosome.contains(sequence.getSequenceName()))
                dictionary.addSequence(sequence);
        }

        return new SAMFileHeader(dictionary);
    }
}
