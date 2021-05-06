package com.hartwig.hmftools.lilac.read;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Collectors;

public class SAMSlicer
{
    private final int mMinMappingQuality;

    public SAMSlicer(int minMappingQuality)
    {
        mMinMappingQuality = minMappingQuality;
    }

    public void slice(final GenomeRegion region, final SamReader samReader, final Consumer<SAMRecord> consumer)
    {
        String string = region.chromosome();
        slice(string, (int) region.start(), (int) region.end(), samReader, consumer);
    }

    public void slice(final BaseRegion region, final SamReader samReader, final Consumer<SAMRecord> consumer)
    {
        slice(region.Chromosome, region.start(), region.end(), samReader, consumer);
    }

    public void slice(final String contig, int start, int end, final SamReader samReader, final Consumer<SAMRecord> consumer)
    {
        SAMRecordIterator sAMRecordIterator = samReader.query(contig, start, end, false);

        while(sAMRecordIterator.hasNext())
        {
            SAMRecord samRecord = sAMRecordIterator.next();

            if(!meetsQualityRequirements(samRecord))
                continue;

            consumer.accept(samRecord);
        }
    }

    public List<SAMRecord> slice(final String contig, int start, int end, final SamReader samReader)
    {
        final List<SAMRecord> records = Lists.newArrayList();

        SAMRecordIterator sAMRecordIterator = samReader.query(contig, start, end, false);

        while(sAMRecordIterator.hasNext())
        {
            SAMRecord samRecord = sAMRecordIterator.next();

            if(!meetsQualityRequirements(samRecord))
                continue;

            records.add(samRecord);
        }

        return records;
    }

    public final List<SAMRecord> queryMates(final SamReader samReader, final List<SAMRecord> records)
    {
        return records.stream().map(x -> samReader.queryMate(x))
                .filter(x -> x != null)
                .filter(x -> meetsQualityRequirements(x))
                .collect(Collectors.toList());
    }

    private final boolean meetsQualityRequirements(SAMRecord record)
    {
        return record.getMappingQuality() >= mMinMappingQuality
            && !record.getReadUnmappedFlag()
            && !record.getDuplicateReadFlag()
            && !record.isSecondaryOrSupplementary();
    }

}
