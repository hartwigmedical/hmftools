package com.hartwig.hmftools.bamtools.compare;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;

import java.util.Comparator;
import java.util.List;

import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;

public class ReadKey
{
    public final String ReadName;
    public final boolean FirstOfPair;

    // supplementary index is 0 for non supplementary, and then rest rank by chr:pos:strand:cigar
    public final byte SupplementaryIndex;

    public ReadKey(final String readName, final boolean firstOfPair, final int supplementaryIndex)
    {
        ReadName = readName;
        FirstOfPair = firstOfPair;
        SupplementaryIndex = (byte)supplementaryIndex;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(this == o)
        {
            return true;
        }
        if(!(o instanceof ReadKey))
        {
            return false;
        }

        final ReadKey readKey = (ReadKey) o;

        if(FirstOfPair != readKey.FirstOfPair)
        {
            return false;
        }
        if(SupplementaryIndex != readKey.SupplementaryIndex)
        {
            return false;
        }
        return ReadName.equals(readKey.ReadName);
    }

    @Override
    public int hashCode()
    {
        int result = ReadName.hashCode();
        result = 31 * result + (FirstOfPair ? 1 : 0);
        result = 31 * result + (int) SupplementaryIndex;
        return result;
    }

    private static final Comparator<SupplementaryReadData> supplementaryAlignmentComparator =
            Comparator.comparing((SupplementaryReadData o) -> o.Chromosome)
            .thenComparingInt(o -> o.Position)
            .thenComparingInt(o -> o.Strand)
            .thenComparing(o -> o.Cigar);

    // create key from read
    public static ReadKey from(SAMRecord read)
    {
        return new ReadKey(read.getReadName(), SamRecordUtils.firstInPair(read), calcSupplementaryIndex(read));
    }

    static int calcSupplementaryIndex(SAMRecord read)
    {
        if(!read.isSecondaryOrSupplementary())
        {
            return 0;
        }

        String saAttribute = read.getStringAttribute(SAMTag.SA.name());

        if(saAttribute == null)
        {
            return 1;
        }

        // before splitting it, we want to see how many supplementary there are, problem is that this split
        // call takes too much memory otherwise. If there is only one supplementary then the index would be 1
        int numSAs = 1;

        // NOTE: we skip the last character as trailing ; is not relevant
        for(int i = 0; i < saAttribute.length() - 1; ++i)
        {
            if(saAttribute.charAt(i) == ';')
            {
                ++numSAs;
            }
        }
        if(numSAs <= 1)
        {
            // first supplementary has index 1
            return 1;
        }
        // split it by ;
        // SA looks like 20,61647163,+,99M52S,0,1;11,70524575,+,95S30M26S,0,0;
        // spec is: SA:Z:(rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+
        // Java split function automatically discard trailing empty elements
        List<SupplementaryReadData> saItems;

        try
        {
            saItems = SupplementaryReadData.extractAlignments(saAttribute);
        }
        catch(IndexOutOfBoundsException | NumberFormatException e)
        {
            BT_LOGGER.warn("unable to parse alignments from SA attribute: {}", saAttribute);
            return 1;
        }

        if(saItems == null || saItems.size() <= 1)
        {
            // only 1 supplementary read, shouldn't happen
            return 1;
        }

        SupplementaryReadData sa = new SupplementaryReadData(read.getReferenceName(), read.getAlignmentStart(),
                read.getReadNegativeStrandFlag() ? '-' : '+', read.getCigarString(), read.getMappingQuality(), 0);

        // first one is the non supplementary read, we skip over that. We only want to rank this read against other
        // supplementaries
        int saIndex = 1;
        for(int i = 1; i < saItems.size(); ++i)
        {
            if(supplementaryAlignmentComparator.compare(sa, saItems.get(i)) > 0)
            {
                ++saIndex;
            }
        }

        return saIndex;
    }
}