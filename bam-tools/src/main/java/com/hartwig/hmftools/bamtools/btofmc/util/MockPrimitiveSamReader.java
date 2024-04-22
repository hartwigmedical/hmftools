package com.hartwig.hmftools.bamtools.btofmc.util;

import static com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome.MT_LENGTH;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.region.BaseRegion;

import org.apache.commons.compress.utils.Lists;
import org.apache.commons.lang3.NotImplementedException;

import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileSpan;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;

public class MockPrimitiveSamReader implements SamReader.PrimitiveSamReader
{
    private static final SAMFileHeader FILE_HEADER;

    static
    {
        SAMSequenceDictionary samDictionaryV37 = new SAMSequenceDictionary();

        RefGenomeCoordinates v37Coords = RefGenomeCoordinates.COORDS_37;

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            samDictionaryV37.addSequence(new SAMSequenceRecord(chromosome.toString(), v37Coords.Lengths.get(chromosome)));
        }

        samDictionaryV37.addSequence(new SAMSequenceRecord(MitochondrialChromosome.MT.toString(), MT_LENGTH));

        FILE_HEADER = new SAMFileHeader();
        FILE_HEADER.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        FILE_HEADER.setSequenceDictionary(samDictionaryV37);
    }

    private final SortedMap<Integer, SortedSet<SAMRecord>> mSamRecordsByRefIdx;
    private final List<SAMRecord> mUnmappedSamRecords;

    public MockPrimitiveSamReader(final Collection<SAMRecord> samRecords)
    {
        mSamRecordsByRefIdx = Maps.newTreeMap();
        mUnmappedSamRecords = Lists.newArrayList();

        populateSamRecords(samRecords);
    }

    private void populateSamRecords(final Collection<SAMRecord> samRecords)
    {
        for(SAMRecord samRecord : samRecords)
        {
            if(samRecord.getReadUnmappedFlag() && samRecord.getMateUnmappedFlag())
            {
                mUnmappedSamRecords.add(samRecord);
                continue;
            }

            int refIdx = samRecord.getReferenceIndex();
            SortedSet<SAMRecord> records = mSamRecordsByRefIdx.get(refIdx);
            if(records == null)
            {
                records = Sets.newTreeSet(Comparator.comparingInt(SAMRecord::getAlignmentStart));
                mSamRecordsByRefIdx.put(refIdx, records);
            }

            records.add(samRecord);
        }
    }

    @Override
    public SamReader.Type type()
    {
        // TODO:
        throw new NotImplementedException("TODO");
    }

    @Override
    public boolean hasIndex()
    {
        throw new NotImplementedException("TODO");
    }

    @Override
    public BAMIndex getIndex()
    {
        throw new NotImplementedException("TODO");
    }

    @Override
    public SAMFileHeader getFileHeader()
    {
        return FILE_HEADER;
    }

    @Override
    public CloseableIterator<SAMRecord> getIterator()
    {
        List<SAMRecord> samRecords = mSamRecordsByRefIdx.values().stream().flatMap(SortedSet::stream).collect(Collectors.toList());
        samRecords.addAll(mUnmappedSamRecords);
        return new IteratorToCloseableIteratorAdaptor<>(samRecords.iterator());
    }

    @Override
    public CloseableIterator<SAMRecord> getIterator(final SAMFileSpan samFileSpan)
    {
        throw new NotImplementedException("TODO");
    }

    @Override
    public SAMFileSpan getFilePointerSpanningReads()
    {
        throw new NotImplementedException("TODO");
    }

    private List<SAMRecord> query(final QueryInterval queryInterval, boolean contained)
    {
        SortedSet<SAMRecord> samRecords = mSamRecordsByRefIdx.get(queryInterval.referenceIndex);
        if(samRecords == null)
        {
            return Lists.newArrayList();
        }

        BaseRegion queryRegion = new BaseRegion(queryInterval.start, queryInterval.end);
        List<SAMRecord> filteredSamRecords = samRecords.stream().filter(samRecord ->
        {
            BaseRegion alignmentRegion = new BaseRegion(samRecord.getAlignmentStart(), samRecord.getAlignmentEnd());
            if(contained)
            {
                return queryRegion.containsRegion(alignmentRegion);
            }

            return queryRegion.overlaps(alignmentRegion);
        }).collect(Collectors.toList());

        return filteredSamRecords;
    }

    @Override
    public CloseableIterator<SAMRecord> query(final QueryInterval[] queryIntervals, boolean contained)
    {
        SortedSet<SAMRecord> samRecords = Sets.newTreeSet(Comparator.comparingInt(SAMRecord::getAlignmentStart));
        for(QueryInterval queryInterval : queryIntervals)
        {
            samRecords.addAll(query(queryInterval, contained));
        }

        return new IteratorToCloseableIteratorAdaptor<>(samRecords.iterator());
    }

    @Override
    public CloseableIterator<SAMRecord> queryAlignmentStart(final String s, final int i)
    {
        throw new NotImplementedException("TODO");
    }

    @Override
    public CloseableIterator<SAMRecord> queryUnmapped()
    {
        return new IteratorToCloseableIteratorAdaptor<>(mUnmappedSamRecords.iterator());
    }

    @Override
    public void close()
    {
    }

    @Override
    public ValidationStringency getValidationStringency()
    {
        return ValidationStringency.SILENT;
    }
}
