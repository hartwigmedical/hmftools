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
    private final boolean mProcessSupplementaries;

    public SAMSlicer(int minMappingQuality, boolean processSupplementaries)
    {
        mMinMappingQuality = minMappingQuality;
        mProcessSupplementaries = processSupplementaries;
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

    public List<SAMRecord> queryMates(final SamReader samReader, final List<SAMRecord> records)
    {
        return records.stream().map(x -> samReader.queryMate(x))
                .filter(x -> x != null)
                .filter(x -> meetsQualityRequirements(x))
                .collect(Collectors.toList());
    }

    private boolean meetsQualityRequirements(SAMRecord record)
    {
        if(record.getMappingQuality() < mMinMappingQuality)
            return false;

        if(record.getReadUnmappedFlag() || record.getDuplicateReadFlag() || record.isSecondaryAlignment())
            return false;

        if(record.getSupplementaryAlignmentFlag() && !mProcessSupplementaries)
            return false;

        return true;
    }

}
