package com.hartwig.hmftools.common.bam;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;

import static htsjdk.samtools.CigarOperator.D;

import java.util.List;

import com.google.common.collect.Lists;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class BamRecordUtils
{
    public static final List<int[]> generateMappedCoords(final Cigar cigar, int posStart)
    {
        final List<int[]> mappedCoords = Lists.newArrayList();

        // first establish whether the read is split across 2 distant regions, and if so which it maps to
        int posOffset = 0;
        boolean continueRegion = false;

        for(CigarElement element : cigar.getCigarElements())
        {
            if(element.getOperator() == CigarOperator.S)
            {
                // nothing to skip
            }
            else if(element.getOperator() == D)
            {
                posOffset += element.getLength();
                continueRegion = true;
            }
            else if(element.getOperator() == CigarOperator.I)
            {
                // nothing to skip
                continueRegion = true;
            }
            else if(element.getOperator() == CigarOperator.N)
            {
                posOffset += element.getLength();
                continueRegion = false;
            }
            else if(element.getOperator() == CigarOperator.M)
            {
                int readStartPos = posStart + posOffset;
                int readEndPos = readStartPos + element.getLength() - 1;

                if(continueRegion && !mappedCoords.isEmpty())
                {
                    int[] lastRegion = mappedCoords.get(mappedCoords.size() - 1);
                    lastRegion[SE_END] = readEndPos;
                }
                else
                {
                    mappedCoords.add(new int[] { readStartPos, readEndPos });
                }

                posOffset += element.getLength();
                continueRegion = false;
            }
        }

        return mappedCoords;
    }

}
