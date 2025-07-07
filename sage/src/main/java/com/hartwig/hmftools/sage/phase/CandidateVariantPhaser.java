package com.hartwig.hmftools.sage.phase;

import static java.lang.Math.log10;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.phase.PhasedVariantGroup.maxPosition;
import static com.hartwig.hmftools.sage.phase.PhasedVariantGroup.minPosition;
import static com.hartwig.hmftools.sage.phase.PhasingUtils.mergeByExtension;
import static com.hartwig.hmftools.sage.phase.PhasingUtils.mergeMatching;
import static com.hartwig.hmftools.sage.phase.PhasingUtils.mergeUninformative;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.perf.PerformanceCounter;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

public class CandidateVariantPhaser implements VariantPhaser
{
    private final PhaseSetCounter mPhaseSetCounter;
    private final boolean mLogData;

    private ChrBaseRegion mRegion; // only for logging

    private final List<PhasedGroupCollection> mPhasedGroupCollections; // order by position of the lowest variant
    private int mNextGroupId;

    private final PerformanceCounter mPerfCounter;

    private static final int INITIAL_MIN_READ_COUNT = 1;
    private static final int FINAL_MIN_READ_COUNT = 2;

    public CandidateVariantPhaser(final PhaseSetCounter phaseSetCounter, boolean logData)
    {
        mPhaseSetCounter = phaseSetCounter;
        mLogData = logData;

        mPhasedGroupCollections = Lists.newArrayList();
        mNextGroupId = 0;

        mPerfCounter = new PerformanceCounter("PhaseReads");
    }

    public List<PhasedGroupCollection> getPhasedCollections() { return mPhasedGroupCollections; }
    public int getPhasingGroupCount() { return mPhasedGroupCollections.stream().mapToInt(x -> x.groupCount()).sum(); }
    public PerformanceCounter getPerfCounter() { return mPerfCounter; }
    public ChrBaseRegion region() { return mRegion; }

    public void clearAll()
    {
        mPhasedGroupCollections.clear();
    }

    @Override
    public void initialise(final ChrBaseRegion region, final String sample)
    {
        mRegion = region;
        mPhasedGroupCollections.clear();
        mNextGroupId = 0;

        mPerfCounter.start();
        mPerfCounter.pause();
    }

    @Override
    public void registeredPhasedVariants(final List<ReadContextCounter> posCounters, final List<ReadContextCounter> negCounters)
    {
        if(posCounters.isEmpty() || posCounters.size() + negCounters.size() < 2)
            return;

        mPerfCounter.resume();
        processPhasedVariants(posCounters, negCounters);
        mPerfCounter.pause();
    }

    private void processPhasedVariants(final List<ReadContextCounter> posCounters, final List<ReadContextCounter> negCounters)
    {
        // walk backwards then forwards from the current location looking for a match
        int posVarMin = minPosition(posCounters, true);
        int posVarMax = maxPosition(posCounters, true);

        int variantMin = !negCounters.isEmpty() ? min(posVarMin, minPosition(negCounters, true)) : posVarMin;
        int variantMax = !negCounters.isEmpty() ? max(posVarMax, maxPosition(negCounters, true)) : posVarMax;

        // assign collections of variants to non-overlapping collections for subsequent comparison efficiency
        int index = 0;
        while(index < mPhasedGroupCollections.size())
        {
            PhasedGroupCollection collection = mPhasedGroupCollections.get(index);
            if(collection.positionsOverlap(variantMin, variantMax))
            {
                if(collection.addPhaseVariants(posVarMin, posVarMax, mNextGroupId, posCounters, negCounters))
                {
                    ++mNextGroupId;

                    // check for a merge with the next collection
                    if(index < mPhasedGroupCollections.size() - 1)
                    {
                        PhasedGroupCollection nextCollection = mPhasedGroupCollections.get(index + 1);
                        if(collection.positionsOverlap(nextCollection.minPosition(), nextCollection.maxPosition()))
                        {
                            collection.merge(nextCollection);
                            mPhasedGroupCollections.remove(index + 1);
                        }
                    }
                }

                return;
            }
            else if(variantMax < collection.minPosition())
            {
                break;
            }

            ++index;
        }

        PhasedGroupCollection collection = new PhasedGroupCollection(mPhasedGroupCollections.size());
        collection.addPhaseVariants(posVarMin, posVarMax, mNextGroupId, posCounters, negCounters);
        ++mNextGroupId;
        mPhasedGroupCollections.add(collection);
    }

    public void signalPhaseReadsEnd() { mPerfCounter.stop(); }

    public void assignLocalPhaseSets(final Set<ReadContextCounter> passingCounters, final Set<ReadContextCounter> validCounters)
    {
        // assign local phase set IDs to all phased variants
        int startCount = mNextGroupId;

        boolean hasGroups = applyInitialFilters(passingCounters, validCounters);

        if(!hasGroups)
            return;

        int startFilteredCount = mPhasedGroupCollections.stream().mapToInt(x -> x.groups().size()).sum();

        if(mLogData)
        {
            List<PhasedVariantGroup> phasedGroups = Lists.newArrayList();
            mPhasedGroupCollections.forEach(x -> phasedGroups.addAll(x.groups()));
            logPhasedReadCounters(phasedGroups, "INITIAL");
        }

        mergeGroups();

        List<PhasedVariantGroup> finalPhasedGroups = Lists.newArrayList();
        mPhasedGroupCollections.forEach(x -> finalPhasedGroups.addAll(x.groups()));

        if(mLogData)
        {
            logPhasedReadCounters(finalPhasedGroups, "FINAL");
        }

        int assignedLps = 0;
        Set<ReadContextCounter> uniqueRCs = mLogData && SG_LOGGER.isTraceEnabled() ? Sets.newHashSet() : null;

        for(PhasedVariantGroup group : finalPhasedGroups)
        {
            if(group.ReadCount + group.AllocatedReadCount < FINAL_MIN_READ_COUNT)
                continue;

            int nextLps = mPhaseSetCounter.getNext();
            group.PositiveReadCounters.forEach(x -> x.addLocalPhaseSet(nextLps, group.ReadCount, group.AllocatedReadCount));
            ++assignedLps;

            if(uniqueRCs != null)
            {
                group.PositiveReadCounters.forEach(x -> uniqueRCs.add(x));
                group.NegativeReadCounters.forEach(x -> uniqueRCs.add(x));
            }
        }

        if(SG_LOGGER.isTraceEnabled())
        {
            SG_LOGGER.trace("region({}) phasing groups(coll={} start={} filtered={}) postMerge({}) assigned({}) rc(pass={} valid={} uniqueRCs={})",
                    mRegion, mPhasedGroupCollections.size(), startCount, startFilteredCount, finalPhasedGroups.size(), assignedLps,
                    passingCounters.size(), validCounters.size(), uniqueRCs != null ? uniqueRCs.size() : 0);
        }

        // mPerfCounters.get(PC_FORM_LPS).stop();
    }

    private boolean applyInitialFilters(
            final Set<ReadContextCounter> passingReadCounters, final Set<ReadContextCounter> validCounters)
    {
        // returns true if a group has variants left after filtering
        if(passingReadCounters.isEmpty())
            return false;

        boolean hasGroups = false;

        for(PhasedGroupCollection collection : mPhasedGroupCollections)
        {
            collection.finalise();

            // require more reads to keep a group where there are high numbers in a collection (ie 2 at 3100, 3 at 31K, 4 at 310K)
            int readCountThreshold = max((int)round(log10(collection.groups().size())) - 2, INITIAL_MIN_READ_COUNT);

            if(readCountThreshold >= 3)
            {
                SG_LOGGER.trace("region({}) phasing group collection({}) sets min read count to {}",
                        mRegion, collection, readCountThreshold);
            }

            int index = 0;
            while(index < collection.groups().size())
            {
                PhasedVariantGroup group = collection.groups().get(index);

                if(group.ReadCount < readCountThreshold)
                {
                    collection.groups().remove(index);
                    continue;
                }

                if(group.PositiveReadCounters.stream().noneMatch(y -> passingReadCounters.contains(y)))
                {
                    collection.groups().remove(index);
                    continue;
                }

                if(group.cullReadCounters(validCounters) && !group.isValid())
                {
                    collection.groups().remove(index);
                    continue;
                }

                ++index;
            }

            hasGroups |= !collection.groups().isEmpty();
        }

        return hasGroups;
    }

    public void mergeGroups()
    {
        // reduce groups to minimum set of information local phase sets by merging:
        // groups with matching positives, non-conflicting negatives and not subsets of other groups
        // extending groups with common overlaps if only one group extends either up or down
        // removing uninformative groups

        for(PhasedGroupCollection collection : mPhasedGroupCollections)
        {
            if(collection.groups().isEmpty())
                continue;

            // then apply merging rules within these overlapping groups
            mergeMatching(collection.groups());

            mergeByExtension(collection.groups());

            mergeUninformative(collection.groups());

            mergeByExtension(collection.groups());
        }
    }

    private void logPhasedReadCounters(final List<PhasedVariantGroup> phasedGroups, final String stage)
    {
        for(PhasedVariantGroup group : phasedGroups)
        {
            StringJoiner posVars = new StringJoiner(";");
            StringJoiner posIds = new StringJoiner(";");
            StringJoiner negVars = new StringJoiner(";");
            StringJoiner negIds = new StringJoiner(";");

            for(ReadContextCounter rc : group.PositiveReadCounters)
            {
                posVars.add(getReadCounterVar(rc));
                posIds.add(String.valueOf(rc.id()));
            }

            for(ReadContextCounter rc : group.NegativeReadCounters)
            {
                negVars.add(getReadCounterVar(rc));
                negIds.add(String.valueOf(rc.id()));
            }

            // Time,Stage,Id,MergedIds,PosIds,NegIds,UniqueReadCount,AllocReadCount,PosVars,NegVars
            SG_LOGGER.debug(String.format("LPS_DATA,%s,%d,%s,%s,%s,%d,%.1f,%s,%s",
                    stage, group.Id, group.mergedGroupIds(), posIds, negIds,
                    group.ReadCount, group.AllocatedReadCount, posVars, negVars));
        }
    }

    private static String getReadCounterVar(final ReadContextCounter rc)
    {
        return String.format("%d_%s>%s", rc.position(), rc.ref(), rc.alt());
    }

}
