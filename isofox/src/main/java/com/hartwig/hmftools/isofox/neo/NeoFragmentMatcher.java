package com.hartwig.hmftools.isofox.neo;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.isofox.neo.NeoFragmentSupport.EXACT_MATCH;
import static com.hartwig.hmftools.isofox.neo.NeoFragmentSupport.MISMATCH;
import static com.hartwig.hmftools.isofox.neo.NeoFragmentSupport.PARTIAL_MATCH;

import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.isofox.common.Read;
import com.hartwig.hmftools.isofox.fusion.ChimericReadGroup;

import com.google.common.collect.Lists;

public class NeoFragmentMatcher
{
    public static final int MIN_BASE_OVERLAP = 10;
    public static final int NOVEL_BASE_OVERLAP = 5;
    private static final int MAX_BASE_MISMATCH = 2;

    public static NeoFragmentSupport findFusionSupport(final NeoEpitopeData neData, int stream, final Read read)
    {
        NeoFragmentSupport support = new NeoFragmentSupport();

        // expect the read to either fully fall within one of the up or down stream ranges, or be soft-clipped
        int junctionSide = neData.Orientations[stream] == ORIENT_FWD ? SE_END : SE_START;

        // reads always go up in position (+ve to -ve orientation)
        final int[] codingBaseRange = neData.getCodingBaseRange(stream);

        // if this is a single-chromosome fusion, the read may extend into the other stream's bases and so support the fusion
        // or it may support an un-fused gene
        boolean readWithinStream = (junctionSide == SE_START && read.alignmentStart() >= codingBaseRange[SE_START])
                || (junctionSide == SE_END && read.alignmentEnd() <= codingBaseRange[SE_END]);

        if(readWithinStream)
        {
            final String neoCodingBases = neData.Source.CodingBases[stream];

            int overlapBases = calcCoordinatesOverlap(read.getMappedRegionCoords(), neData.CodingBaseCoords[stream]);

            if(overlapBases < MIN_BASE_OVERLAP)
                return support;

            int maxStartPos = max(read.alignmentStart(), codingBaseRange[SE_START]);
            int minEndPos = min(read.alignmentEnd(), codingBaseRange[SE_END]);

            int matchLevel = compareCodingBases(read, neoCodingBases, neData.CodingBaseCoords[stream], maxStartPos, minEndPos);

            if(matchLevel == MISMATCH)
                return support;

            // now check if any soft-clipped bases match the bases on the other side of this fusion junction

            // soft-clipped bases from the read which span the fusion junction should match the coding bases on the other stream,
            // after adjusting for strand/orientation

            if(read.isSoftClippedNoRegionMatch(junctionSide))
            {
                final String postJuncCodingBases = neData.getFusionSoftClippedBases(stream);

                int softClipLength = junctionSide == SE_START ? read.leftClipLength() : read.rightClipLength();

                if(softClipLength >= MIN_BASE_OVERLAP / 2)
                {
                    final String readSoftClippedBases = junctionSide == SE_START
                            ?
                            read.readBases().substring(0, softClipLength)
                            : read.readBases().substring(read.baseLength() - softClipLength);

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

            if(neData.Orientations[FS_UP] == ORIENT_FWD)
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
                BaseRegion coordLower = read.getMappedRegionCoords().get(i);
                BaseRegion coordUpper = read.getMappedRegionCoords().get(i + 1);

                if(coordLower.end() == fusionJunction[SE_START] && coordUpper.start() == fusionJunction[SE_END])
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

            int maxStartPos = max(read.alignmentStart(), codingBaseRange[SE_START]);
            int minEndPos = min(read.alignmentEnd(), codingBaseRange[SE_END]);

            int matchLevel = compareCodingBases(read, neoCodingBases, neData.CodingBaseCoords[stream], maxStartPos, minEndPos);

            if(matchLevel == MISMATCH)
                return support;

            if(positionWithin(neData.Positions[FS_UP], read.alignmentStart(), read.alignmentEnd()))
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

    public static int calcCoordinatesOverlap(final List<BaseRegion> baseRegions, final List<int[]> coords2)
    {
        int overlapBases = 0;

        for(BaseRegion region : baseRegions)
        {
            for(int[] coord2 : coords2)
            {
                overlapBases += calcBaseOverlap(region, coord2);
            }
        }

        return overlapBases;
    }

    public static int calcBaseOverlap(final BaseRegion range1, final int[] range2)
    {
        int maxStart = max(range1.start(), range2[SE_START]);
        int minEnd = min(range1.end(), range2[SE_END]);

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
            final Read read, final String neoCodingBases, final List<int[]> neoCoords, int posStart, int posEnd)
    {
        int readBaseIndex = 0;
        String readBases = "";

        if(read.isSoftClippedNoRegionMatch(SE_START))
            readBaseIndex += read.leftClipLength();

        for(BaseRegion mappedCoords : read.getMappedRegionCoords())
        {
            if(posStart > mappedCoords.end())
            {
                readBaseIndex += mappedCoords.end() - mappedCoords.start() + 1;
                continue;
            }

            // will now point at the first base of this next region
            int basePos = 0;
            if(positionWithin(posStart, mappedCoords.start(), mappedCoords.end()))
            {
                readBaseIndex += posStart - mappedCoords.start();
                basePos = posStart;
            }
            else
            {
                basePos = mappedCoords.start();
            }

            for(; basePos <= mappedCoords.end(); ++basePos)
            {
                if(basePos > posEnd || readBaseIndex >= read.baseLength())
                    break;

                readBases += read.readBases().substring(readBaseIndex, readBaseIndex + 1);
                ++readBaseIndex;
            }

            if(basePos > posEnd)
                break;
        }

        int neoBaseIndex = 0;
        String neoBases = "";

        for(int[] neoCoord : neoCoords)
        {
            if(posStart > neoCoord[SE_END])
            {
                neoBaseIndex += neoCoord[SE_END] - neoCoord[SE_START] + 1;
                continue;
            }

            int basePos = 0;
            if(positionWithin(posStart, neoCoord[SE_START], neoCoord[SE_END]))
            {
                neoBaseIndex += posStart - neoCoord[SE_START];
                basePos = posStart;
            }
            else
            {
                basePos = neoCoord[SE_START];
            }

            for(; basePos <= neoCoord[SE_END]; ++basePos)
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
                refBase = neData.getCodingBaseRange(fs)[neData.Orientations[fs] == ORIENT_FWD ? SE_END : SE_START];

            for(Read read : readGroup.reads())
            {
                if(!read.chromosome().equals(chromosome))
                    continue;

                if(read.getMappedRegionCoords().stream().anyMatch(x -> positionWithin(refBase, x.start(), x.end())))
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
