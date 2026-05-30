package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.isofox.common.Read.clippedSide;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.ClippedSide;
import com.hartwig.hmftools.isofox.common.Read;

public class ChimericPosData
{
    public final String Chromosome;
    public final int Position;
    public final String InitialReadId;
    public final String InitialReadBases;

    public int ReadCount;
    public int SuppCount;
    public int DuplicateCount;

    public static final int PosBucketSize = 10;

    public static int positionBucket(int position) { return (int)(floor(position / PosBucketSize) * PosBucketSize); }

    public static String key(final String chromosome, final int posBucket) { return format("%s_%d", chromosome, posBucket); }

    public final List<ChimericRemoteRegion> RemoteRegions;

    public ChimericPosData(final Read read, int positionBucket)
    {
        Chromosome = read.Chromosome;
        Position = positionBucket;
        InitialReadId = read.Id;
        InitialReadBases = read.readBases();
        RemoteRegions = Lists.newArrayList();

        ReadCount = 0;
        SuppCount = 0;
        DuplicateCount = 0;
    }

    public void addReadCounts(final Read read)
    {
        ++ReadCount;

        if(read.isSupplementaryAlignment())
            ++SuppCount;

        if(read.isDuplicate())
            ++DuplicateCount;
    }

    public void addRemoteRegion(final String chromosome, final int regionStart)
    {
        int regionEnd = regionStart + 150; // assumed

        if(Chromosome.equals(chromosome) && positionsOverlap(Position, Position + PosBucketSize, regionStart, regionEnd))
            return;

        ChimericRemoteRegion matchedRegion = RemoteRegions.stream()
                .filter(x -> x.overlaps(chromosome, regionStart, regionEnd)).findFirst().orElse(null);

        if(matchedRegion != null)
        {
            matchedRegion.setStart(min(matchedRegion.start(), regionStart));
            matchedRegion.setEnd(max(matchedRegion.end(), regionEnd));
            ++matchedRegion.Count;
        }
        else
        {
            RemoteRegions.add(new ChimericRemoteRegion(chromosome, regionStart, regionEnd));
        }
    }

    protected static void addChimericPosData(final Map<String,ChimericPosData> chimericPosDataMap, final ChimericReadGroup readGroup)
    {
        // use the read with the longest soft-clip if there is one
        Read primaryRead = null;
        ClippedSide maxClippedSide = null;
        String suppChromosome = "";
        int suppPosition = 0;

        for(Read read : readGroup.reads())
        {
            ClippedSide clippedSide = clippedSide(read);

            if(primaryRead == null)
            {
                primaryRead = read;
                maxClippedSide = clippedSide;
            }
            else if(primaryRead.isSupplementaryAlignment() && !read.isSupplementaryAlignment())
            {
                primaryRead = read;
                maxClippedSide = clippedSide;
            }
            else
            {
                if(clippedSide.Length > maxClippedSide.Length)
                {
                    primaryRead = read;
                    maxClippedSide = clippedSide;
                }
            }

            if(suppChromosome.isEmpty() && read.hasSuppAlignment())
            {
                String[] suppDataItems = read.getSuppAlignment().split(CSV_DELIM, -1);

                if(suppDataItems.length > 2)
                {
                    suppChromosome = suppDataItems[0];
                    suppPosition = Integer.parseInt(suppDataItems[1]);
                }
            }
        }

        int positionBucket = ChimericPosData.positionBucket(primaryRead.PosStart);
        String key = ChimericPosData.key(primaryRead.Chromosome, positionBucket);

        ChimericPosData posData = chimericPosDataMap.get(key);

        if(posData == null)
        {
            posData = new ChimericPosData(primaryRead, positionBucket);
            chimericPosDataMap.put(key, posData);
        }

        posData.addReadCounts(primaryRead);
        posData.addRemoteRegion(primaryRead.mateChromosome(), primaryRead.mateStartPosition());

        if(!suppChromosome.isEmpty())
            posData.addRemoteRegion(suppChromosome, suppPosition);
    }
}
