package com.hartwig.hmftools.isofox.common;

import static com.hartwig.hmftools.isofox.common.CommonUtils.deriveCommonRegions;
import static com.hartwig.hmftools.isofox.common.ReadTranscriptUtils.validTranscriptType;
import static com.hartwig.hmftools.isofox.common.ReadUtils.consensusDuplicateCount;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_INTRON;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.validExonMatch;
import static com.hartwig.hmftools.isofox.common.RegionReadData.NO_EXON;
import static com.hartwig.hmftools.isofox.common.TransMatchType.OTHER_TRANS;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BaseRegion;

public class Fragment
{
    private List<Read> mReads;

    public Fragment(final Read read)
    {
        mReads = List.of(read);
    }

    public Fragment(final Read read1, final Read read2)
    {
        mReads = List.of(read1, read2);
    }

    public List<Read> reads()
    {
        return mReads;
    }

    public String id() { return mReads.get(0).id(); }

    public String chromosome() { return mReads.get(0).chromosome(); }

    // consensus reads from Redux - these are primaries where duplicate reads are also expected
    // to avoid the additional count from the artificially created primary, skip these for any logic which is expression related,
    // (note: only expression uses duplicates)
    public int fragmentCount()
    {
        return 1 + consensusDuplicateCount(mReads.get(0).bamRecord());
    }

    public int minNumLoci()
    {
        return mReads.stream().mapToInt(Read::numLoci).min().getAsInt();
    }

    public List<Read.AltAlignment> altLoci()
    {
        Read minLociRead = mReads.get(0);

        for(Read read : mReads)
        {
            if(read.numLoci() < minLociRead.numLoci())
                minLociRead = read;
        }

        return minLociRead.altLoci();
    }

    public boolean isFullyIntronic()
    {
        return mReads.stream().allMatch(x -> x.getMappedRegions().isEmpty());
    }

    public boolean containsSplit()
    {
        return mReads.stream().anyMatch(Read::containsSplit);
    }

    public int minAlignmentStart()
    {
        return mReads.stream().mapToInt(Read::alignmentStart).min().getAsInt();
    }

    public int maxAlignmentEnd()
    {
        return mReads.stream().mapToInt(Read::alignmentEnd).max().getAsInt();
    }

    public List<BaseRegion> mergedMappings()
    {
        if(mReads.size() == 1)
            return mReads.get(0).getMappedRegionCoords();

        return deriveCommonRegions(mReads.get(0).getMappedRegionCoords(), mReads.get(1).getMappedRegionCoords());
    }

    public List<RegionReadData> uniqueValidRegions()
    {
        List<RegionReadData> validRegions = Lists.newArrayList();

        for(Read read : mReads)
        {
            for(Map.Entry<RegionReadData, RegionMatchType> entry : read.getMappedRegions().entrySet())
            {
                if(validExonMatch(entry.getValue()) && !validRegions.contains(entry.getKey()))
                    validRegions.add(entry.getKey());
            }
        }

        return validRegions;
    }

    public boolean spansMultipleRegions(int transId)
    {
        List<RegionReadData> regions = Lists.newArrayList();

        for(Read read : mReads)
        {
            for(Map.Entry<RegionReadData, RegionMatchType> entry : read.getMappedRegions().entrySet())
            {
                RegionReadData region = entry.getKey();

                if(region.hasTransId(transId) && entry.getValue() != EXON_INTRON && !regions.contains(region))
                {
                    regions.add(region);

                    if(regions.size() > 1)
                        return true;
                }
            }
        }

        return false;
    }

    public List<Integer> validTypeTranscripts()
    {
        List<Integer> transIds = Lists.newArrayList();

        for(int transId : mReads.get(0).getTranscriptClassifications().keySet())
        {
            if(mReads.stream().allMatch(x -> validTranscriptType(x.getTranscriptClassification(transId))))
                transIds.add(transId);
        }

        return transIds;
    }

    public Set<Integer> invalidTranscripts(final List<Integer> validTranscripts)
    {
        Set<Integer> transIds = Sets.newHashSet();

        for(Read read : mReads)
        {
            for(int transId : read.getTranscriptClassifications().keySet())
            {
                if(!validTranscripts.contains(transId))
                    transIds.add(transId);
            }
        }

        return transIds;
    }

    public void setOtherTranscripts(final List<Integer> validTranscripts)
    {
        for(Read read : mReads)
        {
            read.getTranscriptClassifications().entrySet().stream()
                    .filter(x -> validTranscriptType(x.getValue()))
                    .filter(x -> !validTranscripts.contains(x.getKey()))
                    .forEach(x -> x.setValue(OTHER_TRANS));
        }
    }

    public boolean hasTranscriptClassification(int transId, final TransMatchType type)
    {
        return mReads.stream().anyMatch(x -> x.getTranscriptClassification(transId) == type);
    }

    public void trimAdapterBases()
    {
        if(mReads.size() == 2)
            ReadUtils.trimAdapterBases(mReads.get(0), mReads.get(1));
    }

    public boolean readsInDifferentExons(int transId)
    {
        if(mReads.size() < 2)
            return false;

        for(RegionReadData region1 : mReads.get(0).getMappedRegions().keySet())
        {
            int exonRank = region1.getExonRank(transId);

            if(exonRank == NO_EXON)
                continue;

            for(RegionReadData region2 : mReads.get(1).getMappedRegions().keySet())
            {
                if(region2.getExonRank(transId) == exonRank)
                    return false;
            }
        }

        return true;
    }

    public Orientation orientation()
    {
        Read read1 = mReads.get(0);

        if(mReads.size() == 1)
            return read1.orientation();

        Read read2 = mReads.get(1);

        if(read1.orientation() == read2.orientation())
            return null;

        return read1.isFirstOfPair() ? read1.orientation() : read2.orientation();
    }
}