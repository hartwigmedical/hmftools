package com.hartwig.hmftools.lilac.sam;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.function.Consumer;

public class SAMSlicer
{
    private final int mMinMappingQuality;

    public SAMSlicer(int minMappingQuality)
    {
        mMinMappingQuality = minMappingQuality;
    }

    public final void slice(final GenomeRegion region, final SamReader samReader, final Consumer<SAMRecord> consumer)
    {
        String string = region.chromosome();
        slice(string, (int) region.start(), (int) region.end(), samReader, consumer);
    }

    public final void slice(final String contig, int start, int end, final SamReader samReader, final Consumer<SAMRecord> consumer)
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

    public final List<SAMRecord> queryMates(final SamReader samReader, final List<SAMRecord> records)
    {
        return Lists.newArrayList();

        /*
            fun queryMates(samReader: SamReader, records: List<SAMRecord>): List<SAMRecord> {
                return records
                        .mapNotNull { samReader.queryMate(it) }
                        .filter { it.meetsQualityRequirements() }
            }
         */

        /*
        for(SAMRecord record : records)
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

        Iterable $receiver$iv$iv;
        Iterable $receiver$iv = records;
        Iterable iterable = $receiver$iv;
        Collection destination$iv$iv = new ArrayList();
        void $receiver$iv$iv$iv = $receiver$iv$iv;
        Iterator iterator = $receiver$iv$iv$iv.iterator();
        while(iterator.hasNext())
        {
            SAMRecord sAMRecord;
            Object element$iv$iv$iv;
            Object element$iv$iv = element$iv$iv$iv = iterator.next();
            SAMRecord it = (SAMRecord) element$iv$iv;
            boolean bl = false;
            if(samReader.queryMate(it) == null)
            {
                continue;
            }
            SAMRecord it$iv$iv = sAMRecord;
            destination$iv$iv.add(it$iv$iv);
        }
        $receiver$iv = (List) destination$iv$iv;
        $receiver$iv$iv = $receiver$iv;
        destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            SAMRecord it = (SAMRecord) element$iv$iv;
            boolean bl = false;
            if(!meetsQualityRequirements(it))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        return (List) destination$iv$iv;

         */
    }

    private final boolean meetsQualityRequirements(SAMRecord record)
    {
        return record.getMappingQuality() >= mMinMappingQuality && !record.getReadUnmappedFlag()
                && !record.getDuplicateReadFlag() && !record.isSecondaryOrSupplementary();
    }

}
