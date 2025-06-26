package com.hartwig.hmftools.esvee.assembly.phase;

import static java.lang.Math.floor;
import static java.lang.Math.log10;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.assembly.read.BamReader;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.RemoteRegion;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;

import htsjdk.samtools.SAMRecord;

public class RemoteReadExtractor
{
    private final BamReader mBamReader;

    private RemoteRegion mRemoteRegion;
    private final Set<String> mSourceReadIds;
    private final List<Read> mMatchedRemoteReads;

    private int mTotalRemoteReadsSearch;
    private int mTotalRemoteReadsMatched;
    private int mTotalRemoteReadSlices;

    private int mRemoteReadsMatched;
    private int mRemoteReadSlices;

    public RemoteReadExtractor(final BamReader bamReader)
    {
        mBamReader = bamReader;

        mRemoteRegion = null;
        mSourceReadIds = Sets.newHashSet();
        mMatchedRemoteReads = Lists.newArrayList();

        mTotalRemoteReadsSearch = 0;
        mTotalRemoteReadsMatched = 0;
        mTotalRemoteReadSlices = 0;

        mRemoteReadsMatched = 0;
        mRemoteReadSlices = 0;
    }

    public int remoteReadSlices() { return mTotalRemoteReadSlices; }
    public int remoteReadsSearch() { return mTotalRemoteReadsSearch; }
    public int remoteReadsMatched() { return mTotalRemoteReadsMatched; }

    public void resetCounts()
    {
        mTotalRemoteReadsSearch = 0;
        mTotalRemoteReadsMatched = 0;
        mTotalRemoteReadSlices = 0;
    }

    public static List<RemoteRegion> collectCandidateRemoteRegions(
            final JunctionAssembly assembly, final List<JunctionAssembly> phasedAssemblies, boolean checkAssemblyMatches)
    {
        List<RemoteRegion> remoteRegions = Lists.newArrayList();

        // ignore remote regions which overlap an assembly in this phase group with any matching reads
        for(RemoteRegion remoteRegion : assembly.remoteRegions())
        {
            if(remoteRegion.isSuppOnlyRegion())
                continue;

            if(!checkAssemblyMatches)
            {
                remoteRegions.add(remoteRegion);
                continue;
            }

            boolean matchesAssembly = false;

            for(JunctionAssembly otherAssembly : phasedAssemblies)
            {
                if(otherAssembly == assembly)
                    continue;

                if(otherAssembly.discordantOnly()) // allow overlaps with discordant junctions since these don't have precise junctions
                    continue;

                if(!remoteRegion.overlapsAssembly(otherAssembly))
                    continue;

                // otherwise ignore if the overlapping assembly has these remote reads already
                // note the max reads checked to avoid checking all reads in very large groups
                int readsChecked = 0;
                for(String readId : remoteRegion.readIds())
                {
                    if(otherAssembly.support().stream().filter(x -> x.type().isSplitSupport()).anyMatch(x -> x.id().equals(readId)))
                    {
                        matchesAssembly = true;
                        break;
                    }

                    if(otherAssembly.candidateSupport().stream().anyMatch(x -> x.id().equals(readId)))
                    {
                        matchesAssembly = true;
                        break;
                    }

                    ++readsChecked;
                    if(readsChecked > MAX_MATCHED_READ_CHECK)
                        break;
                }

                if(matchesAssembly)
                    break;
            }

            if(!matchesAssembly)
                remoteRegions.add(remoteRegion);
        }

        return remoteRegions;
    }

    private static final int MAX_MATCHED_READ_CHECK = 100;

    private static final int MAX_REGIONS = 100;
    private static final int MAX_TOTAL_READS = 5000;
    private static final int MAX_REGION_READS = 1000;

    public List<Read> extractRemoteRegionReads(final int phaseGroupId, final List<RemoteRegion> remoteRegions, boolean applyThresholds)
    {
        // slices the BAM at the specified locations for mates of reads in one or more assemblies
        // first merges overlapping reads to avoid repeated reads
        RemoteRegion.mergeRegions(remoteRegions);

        Collections.sort(remoteRegions, Comparator.comparingInt(x -> -x.readCount()));

        int minRemoteRegionReadCount = 1;
        int remoteRegionSliceCount = 0;
        int remoteRegionReadCount = 0;

        // reset for this phase group
        mRemoteReadsMatched = 0;
        mRemoteReadSlices = 0;

        int totalRemoteReads = remoteRegions.stream().mapToInt(x -> x.readCount()).sum();

        applyThresholds |= remoteRegions.size() > MAX_REGIONS;

        // impose restrictions on which remote regions are used - excluding those with 1 or very few reads relatively, eg for 1000 then
        // a region will be required to have 2+ reads
        if(applyThresholds && totalRemoteReads >= MAX_REGION_READS)
        {
            minRemoteRegionReadCount = (int)floor(log10(totalRemoteReads));
        }

        List<Read> remoteRegionReads = Lists.newArrayList();

        for(RemoteRegion remoteRegion : remoteRegions)
        {
            if(minRemoteRegionReadCount > 1 && remoteRegion.readCount() < minRemoteRegionReadCount)
                break;

            List<Read> remoteReads = extractRemoteReads(phaseGroupId, remoteRegion, MAX_REGION_READS);

            remoteRegionReads.addAll(remoteReads);

            ++remoteRegionSliceCount;
            remoteRegionReadCount += remoteReads.size();

            if(remoteRegionSliceCount >= MAX_REGIONS || remoteRegionReadCount >= MAX_TOTAL_READS)
                break;
        }

        return remoteRegionReads;
    }

    private List<Read> extractRemoteReads(final int phaseGroupId, final RemoteRegion remoteRegion, int maxReadCount)
    {
        mRemoteRegion = remoteRegion;

        mSourceReadIds.clear();

        if(remoteRegion.readCount() <= maxReadCount)
        {
            mSourceReadIds.addAll(remoteRegion.readIds());
        }
        else
        {
            for(String readId : remoteRegion.readIds())
            {
                mSourceReadIds.add(readId);

                if(mSourceReadIds.size() >= maxReadCount)
                    break;
            }
        }

        mTotalRemoteReadsSearch += mSourceReadIds.size();

        if(mBamReader != null)
        {
            mMatchedRemoteReads.clear();

            SV_LOGGER.trace("remote region({}) slice", mRemoteRegion);

            mBamReader.sliceBam(mRemoteRegion.Chromosome, mRemoteRegion.start(), mRemoteRegion.end(), this::processRecord);

            SV_LOGGER.trace("remote region({}) sourcedReads(matched={} orig={})",
                    mRemoteRegion, mMatchedRemoteReads.size(), mSourceReadIds.size());

            mTotalRemoteReadsMatched += mMatchedRemoteReads.size();
            mRemoteReadsMatched += mMatchedRemoteReads.size();
        }

        ++mTotalRemoteReadSlices;
        ++mRemoteReadSlices;

        if((mRemoteReadSlices % 100) == 0 && mRemoteReadsMatched > 10000)
        {
            SV_LOGGER.debug("pgId({}) remote region read extraction: slices({}) matched({}) last region({}:{}-{})",
                    phaseGroupId, mRemoteReadSlices, mRemoteReadsMatched,
                    mRemoteRegion.Chromosome, mRemoteRegion.start(), mRemoteRegion.end());
        }

        // ignore supplementaries since their bases provide no new assembly sequence information
        return mMatchedRemoteReads.stream().filter(x -> !x.isSupplementary()).collect(Collectors.toList());
    }

    private void processRecord(final SAMRecord record)
    {
        // the read IDs have been trimmed, so has to match on what has been kept
        boolean containedRead = mSourceReadIds.stream().anyMatch(x -> record.getReadName().contains(x));

        if(!containedRead)
            return;

        Read remoteRead = new Read(record);

        if(mBamReader.currentIsReferenceSample())
            remoteRead.markReference();

        mMatchedRemoteReads.add(remoteRead);
    }

    public static void purgeSupplementaryReads(final JunctionAssembly assembly, final List<Read> remoteReads)
    {
        // purge any read which is a supplementary of an existing junction read
        List<SupportRead> junctionSupplementaries = assembly.support()
                .stream().filter(x -> x.supplementaryData() != null).collect(Collectors.toList());

        if(junctionSupplementaries.isEmpty())
            return;

        if(junctionSupplementaries.size() > MAX_REGION_READS || remoteReads.size() > MAX_REGION_READS)
        {
            // remove all remote reads with supp data on the assumption they are linked with a matching assembly read
            int index = 0;
            while(index < remoteReads.size())
            {
                Read read = remoteReads.get(index);

                if(read.supplementaryData() != null)
                    remoteReads.remove(index);
                else
                    ++index;
            }

            return;
        }

        // otherwise look for matches on read ID
        int index = 0;
        while(index < remoteReads.size())
        {
            Read read = remoteReads.get(index);

            SupportRead matchedJunctionRead = junctionSupplementaries
                    .stream().filter(x -> x.matchesFragment(read, false) && x.firstInPair() == read.firstInPair())
                    .findFirst().orElse(null);

            if(matchedJunctionRead != null)
                remoteReads.remove(index);
            else
                ++index;
        }
    }
}
