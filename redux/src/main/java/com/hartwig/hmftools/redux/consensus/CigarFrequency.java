package com.hartwig.hmftools.redux.consensus;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import htsjdk.samtools.SAMRecord;

public class CigarFrequency implements Comparable<CigarFrequency>
{
    public final List<SAMRecord> Reads;
    public final int ClippedLength;

    public CigarFrequency(final SAMRecord read)
    {
        Reads = Lists.newArrayList(read);

        ClippedLength = read.getCigar().getCigarElements().stream()
                .filter(x -> x.getOperator().isClipping()).mapToInt(x -> x.getLength()).sum();
    }

    public String cigar() { return Reads.get(0).getCigarString(); }
    public int frequency() { return Reads.size(); }

    public SAMRecord firstRead() { return Reads.get(0); }
    public String firstReadName() { return Reads.get(0).getReadName(); }

    public void addRead(final SAMRecord read) { Reads.add(read); }

    public void sortReads()
    {
        Collections.sort(Reads, Comparator.comparing(x -> x.getReadName()));
    }

    public String toString()
    {
        return format("cigar(%s) count(%d)", cigar(), frequency());
    }

    @Override
    public int compareTo(final CigarFrequency other)
    {
        if(frequency() != other.frequency())
            return frequency() > other.frequency() ? -1 : 1;

        if(ClippedLength != other.ClippedLength)
            return ClippedLength < other.ClippedLength ? -1 : 1;

        return firstReadName().compareTo(other.firstReadName());
    }

    public static SAMRecord selectTemplateRead(final List<SAMRecord> reads)
    {
        // group read by most frequency CIGARs then clipped length, then mate CIGAR frequency, then read ID alphabetically
        Map<String,CigarFrequency> cigarFrequencies = buildCigarFrequencies(reads, false);

        if(cigarFrequencies.size() == 1)
        {
            return selectTemplateRead(cigarFrequencies.values().iterator().next());
        }

        List<CigarFrequency> frequencies = cigarFrequencies.values().stream().collect(Collectors.toList());
        Collections.sort(frequencies);

        return selectTemplateRead(frequencies.get(0));
    }

    public static SAMRecord selectTemplateRead(final CigarFrequency cigarFrequency)
    {
        // select read by most common mate cigar
        Map<String,CigarFrequency> mateCigarFrequencies = buildCigarFrequencies(cigarFrequency.Reads, true);

        if(mateCigarFrequencies.size() == 1)
        {
            return mateCigarFrequencies.values().iterator().next().firstRead();
        }

        List<CigarFrequency> mateFrequencies = mateCigarFrequencies.values().stream().collect(Collectors.toList());
        Collections.sort(mateFrequencies);

        return mateFrequencies.get(0).firstRead();
    }

    private static Map<String,CigarFrequency> buildCigarFrequencies(final List<SAMRecord> reads, boolean useMateCigar)
    {
        Map<String,CigarFrequency> cigarFrequencies = Maps.newHashMap();

        for(SAMRecord read : reads)
        {
            final String cigarStr = useMateCigar ? read.getStringAttribute(MATE_CIGAR_ATTRIBUTE) : read.getCigarString();

            CigarFrequency frequency = cigarFrequencies.get(cigarStr);

            if(frequency == null)
                cigarFrequencies.put(cigarStr, new CigarFrequency(read));
            else
                frequency.addRead(read);
        }

        cigarFrequencies.values().forEach(x -> x.sortReads());

        return cigarFrequencies;
    }

}
