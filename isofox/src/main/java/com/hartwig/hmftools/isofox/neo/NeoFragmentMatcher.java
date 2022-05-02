package com.hartwig.hmftools.isofox.neo;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.VAR_INFO_DELIM;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.isofox.neo.NeoFragmentSupport.EXACT_MATCH;
import static com.hartwig.hmftools.isofox.neo.NeoFragmentSupport.MISMATCH;
import static com.hartwig.hmftools.isofox.neo.NeoFragmentSupport.PARTIAL_MATCH;

import java.util.List;

import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.fusion.ChimericReadGroup;

import org.apache.commons.compress.utils.Lists;

public class NeoFragmentMatcher
{
    public static final int MIN_BASE_OVERLAP = 10;
    public static final int NOVEL_BASE_OVERLAP = 5;
    private static final int MAX_BASE_MISMATCH = 2;

    public static NeoFragmentSupport findPointMutationSupport(final NeoEpitopeData neData, final ReadRecord read)
    {
        NeoFragmentSupport support = new NeoFragmentSupport();

        // reads always go up in position (+ve to -ve orientation)
        final int[] codingBaseRange = neData.getCodingBaseRange(FS_UP);
        final String neoCodingBases = neData.getFullCodingBases(FS_UP);

        // reads with soft-clipped bases to another gene will match the coding bases starting midway through the coding bases
        // the read may cross any part of of the coding base region, which will cover all neo-epitope bases (up, down & novel)
        int overlapBases = calcCoordinatesOverlap(read.getMappedRegionCoords(), neData.CodingBaseCoords[FS_UP]);

        if(overlapBases < MIN_BASE_OVERLAP)
            return support;

        int mutationPosition = neData.Positions[FS_UP];
        int maxStartPos = max(read.PosStart, codingBaseRange[SE_START]);
        int minEndPos = min(read.PosEnd, codingBaseRange[SE_END]);

        int matchLevel = compareCodingBases(read, neoCodingBases, neData.CodingBaseCoords[FS_UP], maxStartPos, minEndPos);

        if(matchLevel == MISMATCH)
            return support;

        final String[] varData = neData.Source.VariantInfo.split(VAR_INFO_DELIM);
        int refLength = varData[2].length();
        int altLength = varData[3].length();
        int baseDiff = altLength - refLength;
        int novelLength = 0;

        if(baseDiff == 0)
        {
            novelLength = refLength + NOVEL_BASE_OVERLAP;
        }
        else if(baseDiff > 0)
        {
            // insert case - limit to 10 inserted bases
            novelLength = min(altLength, 10) + NOVEL_BASE_OVERLAP;
        }
        else
        {
            // delete
            novelLength = refLength + NOVEL_BASE_OVERLAP;
        }

        List<int[]> novelRanges = Lists.newArrayList();
        novelRanges.add(new int[] { mutationPosition, mutationPosition });
        expandRange(novelRanges, mutationPosition, neData.CodingBaseCoords[FS_UP], NOVEL_BASE_OVERLAP, false);
        expandRange(novelRanges, mutationPosition, neData.CodingBaseCoords[FS_UP], novelLength, true);

        int novelOverlap = calcCoordinatesOverlap(read.getMappedRegionCoords(), novelRanges);

        if(novelOverlap >= MIN_BASE_OVERLAP)
        {
            ++support.NovelFragments[matchLevel];
        }
        else
        {
            int novelRangeStart = novelRanges.get(0)[SE_START];
            int novelRangeEnd = novelRanges.get(novelRanges.size() - 1)[SE_END];

            if((read.PosEnd < novelRangeEnd && neData.posStrand(FS_UP)) || (read.PosStart > novelRangeStart && !neData.posStrand(FS_UP)))
                ++support.UpFragments[matchLevel];
            else
                ++support.DownFragments[matchLevel];
        }

        return support;
    }

    public static NeoFragmentSupport findFusionSupport(final NeoEpitopeData neData, int stream, final ReadRecord read)
    {
        NeoFragmentSupport support = new NeoFragmentSupport();

        // expect the read to either fully fall within one of the up or down stream ranges, or be soft-clipped
        int junctionSide = neData.Orientations[stream] == POS_ORIENT ? SE_END : SE_START;

        // reads always go up in position (+ve to -ve orientation)
        final int[] codingBaseRange = neData.getCodingBaseRange(stream);

        // if this is a single-chromosome fusion, the read may extend into the other stream's bases and so support the fusion
        // or it may support an un-fused gene
        boolean readWithinStream = (junctionSide == SE_START && read.PosStart >= codingBaseRange[SE_START])
                || (junctionSide == SE_END && read.PosEnd <= codingBaseRange[SE_END]);

        if(readWithinStream)
        {
            final String neoCodingBases = neData.Source.CodingBases[stream];

            int overlapBases = calcCoordinatesOverlap(read.getMappedRegionCoords(), neData.CodingBaseCoords[stream]);

            if(overlapBases < MIN_BASE_OVERLAP)
                return support;

            int maxStartPos = max(read.PosStart, codingBaseRange[SE_START]);
            int minEndPos = min(read.PosEnd, codingBaseRange[SE_END]);

            int matchLevel = compareCodingBases(read, neoCodingBases, neData.CodingBaseCoords[stream], maxStartPos, minEndPos);

            if(matchLevel == MISMATCH)
                return support;

            // now check if any soft-clipped bases match the bases on the other side of this fusion junction

            // soft-clipped bases from the read which span the fusion junction should match the coding bases on the other stream,
            // after adjusting for strand/orientation

            if(read.isSoftClipped(junctionSide))
            {
                final String postJuncCodingBases = neData.getFusionSoftClippedBases(stream);
                int softClipLength = junctionSide == SE_START ?
                        read.Cigar.getFirstCigarElement().getLength() : read.Cigar.getLastCigarElement().getLength();

                if(softClipLength >= MIN_BASE_OVERLAP / 2)
                {
                    final String readSoftClippedBases = junctionSide == SE_START
                            ?
                            read.ReadBases.substring(0, softClipLength)
                            : read.ReadBases.substring(read.ReadBases.length() - softClipLength);

                    matchLevel = calcBaseMatch(postJuncCodingBases, readSoftClippedBases);

                    if(matchLevel == MISMATCH)
                        return support;

                    ++support.NovelFragments[matchLevel];
                }
            }
            else
            {
                if(stream == FS_UP)
                    ++support.UpFragments[matchLevel];
                else
                    ++support.DownFragments[matchLevel];
            }
        }
        else
        {
            if(!neData.isDeletionFusion()) // needs to be soft-clipped if spanning chromosomes
                return support;

            final String neoCodingBases = neData.getFullCodingBases(stream);

            // combine the coords
            final List<int[]> codingBaseCoords = Lists.newArrayList();
            final int[] fusionJunction = new int[SE_PAIR];

            if(neData.Orientations[FS_UP] == POS_ORIENT)
            {
                fusionJunction[SE_START] = neData.Source.CodingBasePositions[FS_UP][SE_END];
                fusionJunction[SE_END] = neData.Source.CodingBasePositions[FS_DOWN][SE_START];
                codingBaseCoords.addAll(neData.CodingBaseCoords[FS_UP]);
                codingBaseCoords.addAll(neData.CodingBaseCoords[FS_DOWN]);
            }
            else
            {
                fusionJunction[SE_START] = neData.Source.CodingBasePositions[FS_DOWN][SE_END];
                fusionJunction[SE_END] = neData.Source.CodingBasePositions[FS_UP][SE_START];
                codingBaseCoords.addAll(neData.CodingBaseCoords[FS_DOWN]);
                codingBaseCoords.addAll(neData.CodingBaseCoords[FS_UP]);
            }

            // the read must have an N-split matching the fusion junction
            boolean supportsSplit = false;

            for(int i = 0; i < read.getMappedRegionCoords().size() - 1; ++i)
            {
                final int[] coordLower = read.getMappedRegionCoords().get(i);
                final int[] coordUpper = read.getMappedRegionCoords().get(i + 1);

                if(coordLower[SE_END] == fusionJunction[SE_START] && coordUpper[SE_START] == fusionJunction[SE_END])
                {
                    supportsSplit = true;
                    break;
                }
            }

            if(!supportsSplit)
                return support;

            int overlapBases = calcCoordinatesOverlap(read.getMappedRegionCoords(), neData.CodingBaseCoords[stream]);

            if(overlapBases < MIN_BASE_OVERLAP)
                return support;

            int maxStartPos = max(read.PosStart, codingBaseRange[SE_START]);
            int minEndPos = min(read.PosEnd, codingBaseRange[SE_END]);

            int matchLevel = compareCodingBases(read, neoCodingBases, neData.CodingBaseCoords[stream], maxStartPos, minEndPos);

            if(matchLevel == MISMATCH)
                return support;

            if(positionWithin(neData.Positions[FS_UP], read.PosStart, read.PosEnd))
            {
                ++support.NovelFragments[matchLevel];
            }
            else
            {
                if(stream == FS_UP)
                    ++support.UpFragments[matchLevel];
                else
                    ++support.DownFragments[matchLevel];
            }
        }

        return support;
    }

    public static int calcCoordinatesOverlap(final List<int[]> coords1, final List<int[]> coords2)
    {
        int overlapBases = 0;

        for(int[] coord1 : coords1)
        {
            for(int[] coord2 : coords2)
            {
                overlapBases += calcBaseOverlap(coord1, coord2);
            }
        }

        return overlapBases;
    }

    public static int calcBaseOverlap(final int[] range1, final int[] range2)
    {
        int maxStart = max(range1[SE_START], range2[SE_START]);
        int minEnd = min(range1[SE_END], range2[SE_END]);

        return maxStart <= minEnd ? minEnd - maxStart + 1 : 0;
    }

    public static void expandRange(final List<int[]> ranges, int position, final List<int[]> coordsList, int shiftCount, boolean shiftUp)
    {
        int coordIndex = 0;

        for(; coordIndex < coordsList.size(); ++coordIndex)
        {
            final int[] coords = coordsList.get(coordIndex);

            if(!positionWithin(position, coords[SE_START], coords[SE_END]))
                continue;

            break;
        }

        if(coordIndex >= coordsList.size())
            return;

        int[] currentRange = !shiftUp ? ranges.get(0) : ranges.get(ranges.size() - 1);

        int currentPos = position;
        int shiftedBases = 0;
        int[] currentCoords = coordsList.get(coordIndex);

        while(shiftedBases < shiftCount)
        {
            if(!shiftUp)
            {
                if(currentPos <= currentCoords[SE_START])
                {
                    --coordIndex;

                    if(coordIndex < 0)
                        return;

                    currentCoords = coordsList.get(coordIndex);
                    currentPos = currentCoords[SE_END];

                    currentRange = new int[] { currentPos, currentPos };
                    ranges.add(0, currentRange);
                }
                else
                {
                    --currentPos;
                }
            }
            else
            {
                if(currentPos >= currentCoords[SE_END])
                {
                    ++coordIndex;

                    if(coordIndex >= coordsList.size())
                        return;

                    currentCoords = coordsList.get(coordIndex);
                    currentPos = currentCoords[SE_START];

                    currentRange = new int[] { currentPos, currentPos };
                    ranges.add(currentRange);

                }
                else
                {
                    ++currentPos;
                }
            }

            ++shiftedBases;
        }

        if(shiftUp)
            currentRange[SE_END] = currentPos;
        else
            currentRange[SE_START] = currentPos;
    }

    public static int compareCodingBases(
            final ReadRecord read, final String neoCodingBases, final List<int[]> neoCoords,
            int posStart, int posEnd)
    {
        int readBaseIndex = 0;
        String readBases = "";

        if(read.isSoftClipped(SE_START))
            readBaseIndex += read.Cigar.getFirstCigarElement().getLength();

        for(int[] mappedCoords : read.getMappedRegionCoords())
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

        return calcBaseMatch(readBases, neoBases);
    }

    private static int calcBaseMatch(final String bases1, final String bases2)
    {
        if(bases1.equals(bases2) || bases1.contains(bases2) || bases2.contains(bases1))
        {
            return EXACT_MATCH;
        }
        else
        {
            // ISF_LOGGER.debug("inexact match");
            int mismatches = 0;
            int overlapBases = min(bases1.length(), bases2.length());

            for(int base = 0; base < overlapBases; ++base)
            {
                char readBase = bases1.charAt(base);
                char neoBase = bases2.charAt(base);

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

    public static void checkBaseCoverage(final NeoEpitopeData neData, final ChimericReadGroup readGroup)
    {
        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            if(neData.isPointMutation() && fs == FS_DOWN)
            {
                // report as the same
                neData.getFragmentSupport().RefBaseDepth[FS_DOWN] = neData.getFragmentSupport().RefBaseDepth[FS_UP];
                break;
            }

            boolean coversBase = false;
            final String chromosome = neData.Chromosomes[fs];

            int refBase;

            if(neData.isPointMutation())
                refBase = neData.Positions[fs];
            else
                refBase = neData.getCodingBaseRange(fs)[neData.Orientations[fs] == POS_ORIENT ? SE_END : SE_START];

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
                ++neData.getFragmentSupport().RefBaseDepth[fs];
        }
    }


}
