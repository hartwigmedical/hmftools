package com.hartwig.hmftools.isofox.neo;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.neo.NeoFragmentSupport.EXACT_MATCH;
import static com.hartwig.hmftools.isofox.neo.NeoFragmentSupport.MISMATCH;
import static com.hartwig.hmftools.isofox.neo.NeoFragmentSupport.PARTIAL_MATCH;

import java.util.List;

import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.fusion.ReadGroup;

public class NeoFragmentMatcher
{
    public static final int MIN_BASE_OVERLAP = 10;
    private static final int MAX_BASE_MISMATCH = 2;

    public static NeoFragmentSupport getNeoEpitopeSupport(final NeoEpitopeData neData, int stream, final ReadRecord read)
    {
        NeoFragmentSupport support = new NeoFragmentSupport();

        // reads always go up in position (+ve to -ve orientation)
        final int[] codingBaseRange = neData.getCodingBaseRange(stream);
        final String neoCodingBases = neData.getFullCodingBases(stream);

        // reads with soft-clipped bases to another gene will match the coding bases starting midway through the coding bases
        if(neData.isPointMutation())
        {
            // the read may cross any part of of the coding base region, which will cover all neo-epitope bases (up, down & novel)
            int maxStartPos = max(read.PosStart, codingBaseRange[SE_START]);
            int minEndPos = min(read.PosEnd, codingBaseRange[SE_END]);

            if(minEndPos - maxStartPos < MIN_BASE_OVERLAP - 1)
                return support;

            int matchLevel = compareCodingBases(read, neoCodingBases, neData.CodingBaseCoords[stream], maxStartPos, minEndPos);

            if(matchLevel == MISMATCH)
                return support;

            int mutationPosition = neData.Positions[FS_UP];

            int[] novelRange = new int[] {
                    mutationPosition - MIN_BASE_OVERLAP / 2,
                    mutationPosition + neData.Source.NovelAA.length() * 3 + MIN_BASE_OVERLAP / 2 };

            int novelOverlap = read.getMappedRegionCoords().stream().map(x -> calcBaseOverlap(x, novelRange)).mapToInt(x -> x).max().orElse(0);

            if(novelOverlap >= MIN_BASE_OVERLAP)
            {
                ++support.NovelFragments[matchLevel];
            }
            else if((read.PosEnd < novelRange[SE_END] && neData.posStrand(FS_UP))
            || (read.PosStart > novelRange[SE_START] && !neData.posStrand(FS_UP)))
            {
                ++support.UpFragments[matchLevel];
            }
            else
            {
                ++support.DownFragments[matchLevel];
            }
        }
        else
        {

        }

        // expect the read to either fully fall within one of the up or down stream ranges, or be soft-clippedf

        return support;
    }

    private static int calcBaseOverlap(final int[] range1, final int[] range2)
    {
        int maxStart = max(range1[SE_START], range2[SE_START]);
        int minEnd = min(range1[SE_END], range2[SE_END]);
        return minEnd - maxStart;
    }

    public static int compareCodingBases(
            final ReadRecord read, final String neoCodingBases, final List<int[]> neoCoords,
            int posStart, int posEnd)
    {
        int readBaseIndex = 0;
        String readBases = "";

        if(read.isSoftClipped(SE_START))
            readBaseIndex += read.Cigar.getFirstCigarElement().getLength();

        for(int[] mappedCoords : read.getMappedRegionCoords(false))
        {
            if(posStart > mappedCoords[SE_END])
            {
                readBaseIndex += mappedCoords[SE_END] - mappedCoords[SE_START] + 1;
                continue;
            }

            // will now point at the first base of this next region
            int basePos = 0;
            if(positionWithin(posStart, mappedCoords[SE_START], mappedCoords[SE_END]))
            {
                readBaseIndex += posStart - mappedCoords[SE_START];
                basePos = posStart;
            }
            else
            {
                basePos = mappedCoords[SE_START];
            }

            for(; basePos <= mappedCoords[SE_END]; ++basePos)
            {
                if(basePos > posEnd || readBaseIndex >= read.ReadBases.length())
                    break;

                readBases += read.ReadBases.substring(readBaseIndex, readBaseIndex + 1);
                ++readBaseIndex;
            }

            if(basePos > posEnd)
                break;
        }

        int neoBaseIndex = 0;
        String neoBases = "";

        for(int[] mappedCoords : neoCoords)
        {
            if(posStart > mappedCoords[SE_END])
            {
                neoBaseIndex += mappedCoords[SE_END] - mappedCoords[SE_START] + 1;
                continue;
            }

            int basePos = 0;
            if(positionWithin(posStart, mappedCoords[SE_START], mappedCoords[SE_END]))
            {
                neoBaseIndex += posStart - mappedCoords[SE_START];
                basePos = posStart;
            }
            else
            {
                basePos = mappedCoords[SE_START];
            }

            for(; basePos <= mappedCoords[SE_END]; ++basePos)
            {
                if(basePos > posEnd || neoBaseIndex >= neoCodingBases.length())
                    break;

                neoBases += neoCodingBases.substring(neoBaseIndex, neoBaseIndex + 1);
                ++neoBaseIndex;
            }

            if(basePos > posEnd)
                break;
        }

        if(readBases.equals(neoBases))
        {
            return EXACT_MATCH;
        }
        else
        {
            // ISF_LOGGER.debug("inexact match");
            int mismatches = 0;
            int overlapBases = min(readBases.length(), neoBases.length());

            for(int base = 0; base < overlapBases; ++base)
            {
                char readBase = readBases.charAt(base);
                char neoBase = neoBases.charAt(base);

                if(readBase != neoBase)
                {
                    ++mismatches;

                    if(mismatches > MAX_BASE_MISMATCH)
                        return MISMATCH;
                }
            }

            return PARTIAL_MATCH;
        }
    }

    public static void checkBaseCoverage(final NeoEpitopeData neData, final ReadGroup readGroup)
    {
        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            if(neData.isPointMutation() && fs == FS_DOWN)
                break;

            boolean coversBase = false;
            final String chromosome = neData.Chromosomes[fs];
            int refBase = neData.Positions[fs];

            for(ReadRecord read : readGroup.Reads)
            {
                if(!read.Chromosome.equals(chromosome))
                    continue;

                if(read.getMappedRegionCoords().stream().anyMatch(x -> positionWithin(refBase, x[SE_START], x[SE_END])))
                {
                    coversBase = true;
                    break;
                }
            }

            if(coversBase)
                ++neData.getFragmentSupport().RefBaseSupport[fs];
        }
    }


}
