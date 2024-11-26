package com.hartwig.hmftools.redux.consensus;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.redux.common.FragmentCoords;

import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;

public class TemplateReadData
{
    public final String ReadId;
    public final String Chromosome;
    public final int AlignmentStart;
    public final String Cigar;
    public final int Flags;

    public TemplateReadData(final String readId, final String chromosome, final int alignmentStart, final String cigar, final int flags)
    {
        ReadId = readId;
        Chromosome = chromosome;
        AlignmentStart = alignmentStart;
        Cigar = cigar;
        Flags = flags;
    }

    public boolean firstInPair() { return isFlagSet(SAMFlag.FIRST_OF_PAIR); }
    public boolean readNegativeStrandFlag() { return isFlagSet(SAMFlag.READ_REVERSE_STRAND); }
    public boolean mateNegativeStrandFlag() { return isFlagSet(SAMFlag.MATE_REVERSE_STRAND); }

    private boolean isFlagSet(final SAMFlag flag) { return (Flags & flag.intValue()) != 0; }

    public static TemplateReadData fromRead(final SAMRecord read)
    {
        return new TemplateReadData(
                read.getReadName(), read.getReferenceName(), read.getAlignmentStart(), read.getCigarString(), read.getFlags());
    }

    public static SAMRecord selectTemplateRead(final List<SAMRecord> reads)
    {
        // establish if the read, its mate or the supplementary data should be used to establish the template data
        FragmentCoords fragmentCoords = FragmentCoords.fromRead(reads.get(0));

        int primaryCount = reads.stream().mapToInt(x -> x.getSupplementaryAlignmentFlag() ? 0 : 1).sum();
        boolean arePrimaries;

        if(primaryCount == reads.size() / 2)
        {
            // select based on the first read's primary status
            Collections.sort(reads, Comparator.comparing(x -> x.getReadName()));
            arePrimaries = !reads.get(0).getSupplementaryAlignmentFlag();
        }
        else
        {
            arePrimaries = primaryCount > reads.size() / 2;
        }

        boolean isLowerRead = fragmentCoords.ReadIsLower;

        List<ReadCigarInfo> readCigarInfos = reads.stream()
                .map(x -> new ReadCigarInfo(x, isLowerRead, arePrimaries)).collect(Collectors.toList());

        String topReadCigar = findTopFrequencyCigar(readCigarInfos, false);

        List<ReadCigarInfo> filteredCigarInfos = readCigarInfos.stream().filter(x -> x.LowerCigar.equals(topReadCigar)).collect(Collectors.toList());

        String topMateCigar = findTopFrequencyCigar(readCigarInfos, true);

        filteredCigarInfos = readCigarInfos.stream().filter(x -> x.MateCigar.equals(topMateCigar)).collect(Collectors.toList());

        if(filteredCigarInfos.size() > 1)
        {
            Collections.sort(filteredCigarInfos, Comparator.comparing(x -> x.Read.getReadName()));
        }

        return filteredCigarInfos.get(0).Read;
    }

    private static String findTopFrequencyCigar(final List<ReadCigarInfo> readCigarInfos, boolean useMateCigar)
    {
        Map<String,Integer> cigarFrequencies = Maps.newHashMap();

        for(ReadCigarInfo readInfo : readCigarInfos)
        {
            String cigarStr = useMateCigar ? readInfo.MateCigar : readInfo.LowerCigar;

            Integer frequency = cigarFrequencies.get(cigarStr);
            cigarFrequencies.put(cigarStr, frequency != null ? frequency + 1 : 1);
        }

        int maxFrequency = cigarFrequencies.values().stream().mapToInt(x -> x.intValue()).max().orElse(0);

        List<String> maxFrequencyCigars = cigarFrequencies.entrySet().stream()
                .filter(x -> x.getValue() == maxFrequency).map(x -> x.getKey()).collect(Collectors.toList());

        Collections.sort(maxFrequencyCigars);
        return maxFrequencyCigars.get(0);
    }
}
