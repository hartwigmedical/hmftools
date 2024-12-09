package com.hartwig.hmftools.redux.consensus;

import static com.hartwig.hmftools.common.bam.CigarUtils.calcCigarAlignedLength;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.redux.common.FragmentCoords;

import htsjdk.samtools.SAMRecord;

public final class TemplateReads
{
    public static SAMRecord selectTemplateRead(final List<SAMRecord> reads, final FragmentCoords fragmentCoords)
    {
        // establish if the read, its mate or the supplementary data should be used to establish the template data
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

        String topMateCigar = findTopFrequencyCigar(filteredCigarInfos, true);

        filteredCigarInfos = filteredCigarInfos.stream().filter(x -> x.MateCigar.equals(topMateCigar)).collect(Collectors.toList());

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
            cigarFrequencies.merge(cigarStr, 1, Integer::sum);
        }

        if(cigarFrequencies.size() == 1)
            return useMateCigar ? readCigarInfos.get(0).MateCigar : readCigarInfos.get(0).LowerCigar;

        int maxFrequency = cigarFrequencies.values().stream().mapToInt(x -> x.intValue()).max().orElse(0);

        List<String> maxFrequencyCigars = cigarFrequencies.entrySet().stream()
                .filter(x -> x.getValue() == maxFrequency).map(x -> x.getKey()).collect(Collectors.toList());

        if(maxFrequencyCigars.size() == 1)
            return maxFrequencyCigars.get(0);

        Collections.sort(maxFrequencyCigars);

        // take the cigar with the most aligned bases
        String maxAlignerCigar = "";
        int maxAligedBases = 0;

        for(String cigarStr : maxFrequencyCigars)
        {
            int alignedLength = cigarStr.equals(NO_CIGAR) ? 1 : calcCigarAlignedLength(cigarStr);

            if(alignedLength > maxAligedBases)
            {
                maxAlignerCigar = cigarStr;
                maxAligedBases = alignedLength;
            }
        }

        return maxAlignerCigar;
    }
}
