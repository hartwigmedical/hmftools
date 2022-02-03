package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;

public class VariantPhaser
{
    private final PhaseSetCounter mPhaseSetCounter;

    private boolean mEnabled;

    private final Map<String,PhasedReadCounters> mPhasedReadCounters;

    public static final int MIN_READ_COUNT = 1;

    public VariantPhaser(final PhaseSetCounter phaseSetCounter)
    {
        mPhaseSetCounter = phaseSetCounter;
        mEnabled = false;

        mPhasedReadCounters = Maps.newHashMap();
    }

    public void setEnabled(boolean toggle) { mEnabled = toggle; }
    public boolean enabled() { return mEnabled; }

    public void reset()
    {
        mPhasedReadCounters.clear();
    }

    public void registeredPhasedVariants(final List<ReadContextCounter> posCounters, final List<ReadContextCounter> negCounters)
    {
        if(!mEnabled)
            return;

        if(posCounters.isEmpty() || posCounters.size() + negCounters.size() < 2)
            return;

        String phaseRcId = formRcId(posCounters, negCounters);

        PhasedReadCounters phasedRc = mPhasedReadCounters.get(phaseRcId);

        if(phasedRc != null)
        {
            phasedRc.ReadCount++;
        }
        else
        {
            int id = mPhasedReadCounters.size();
            mPhasedReadCounters.put(phaseRcId, new PhasedReadCounters(id, phaseRcId, posCounters, negCounters));
        }
    }

    private void mergeGroups()
    {
        /*
        For each additional row update the set of  candidate LPS using the following rules:  (this will make a maximally collapsed set of candidates)
            If it is a subset of an existing LPS then do nothing
            If it agrees with but extends that LPS then extend the LPS
            If it does not match any overlapping LPS then make a new LPS
        For each LPS, count uniquely matching reads and assigned matching reads (pro rata)
        Filter for LPS support > X
        Remove any uninformative LPS (ie LPS has 1+ variant and is the only LPS that includes that variant).
        */

        List<PhasedReadCounters> removedGroups = Lists.newArrayList();

        for(PhasedReadCounters phasedRc : mPhasedReadCounters.values())
        {
            if(removedGroups.contains(phasedRc))
                continue;

            // check for super-set groups
            List<PhasedReadCounters> superGroups = mPhasedReadCounters.values().stream()
                    .filter(x -> x != phasedRc)
                    .filter(x -> !removedGroups.contains(x))
                    .filter(x -> phasedRc.isSubsetOf(x))
                    .collect(Collectors.toList());

            if(superGroups.isEmpty())
                continue;

            // where the super groups differ, allocate by pro-rata and continue
            if(superGroups.stream().anyMatch(x -> x.PositiveReadCounters.size() > phasedRc.PositiveReadCounters.size()))
            {
                double totalReads = superGroups.stream().mapToInt(x -> x.ReadCount).sum();

                for(PhasedReadCounters superGroup : superGroups)
                {
                    double allocFraction = superGroup.ReadCount / totalReads;
                    superGroup.AllocatedReadCount += allocFraction * phasedRc.ReadCount;
                    superGroup.AllocatedReadCount += allocFraction * phasedRc.AllocatedReadCount;
                }

                // remove subset group
                removedGroups.add(phasedRc);
            }
            else
            {
                // otherwise collapse the groups into this one
                for(PhasedReadCounters superGroup : superGroups)
                {
                    phasedRc.ReadCount += superGroup.ReadCount;
                    phasedRc.AllocatedReadCount += superGroup.AllocatedReadCount;
                }

                removedGroups.addAll(superGroups);
            }
        }

        removedGroups.forEach(x -> mPhasedReadCounters.remove(x.RcId));
        removedGroups.clear();

        // remove any small groups or those with a single variant not present as a +ve in another group
        for(PhasedReadCounters phasedRc : mPhasedReadCounters.values())
        {
            if(phasedRc.PositiveReadCounters.size() != 1)
                continue;

            if(mPhasedReadCounters.values().stream()
                    .filter(x -> x != phasedRc)
                    .filter(x -> !removedGroups.contains(x))
                    .noneMatch(x -> x.PositiveReadCounters.contains(phasedRc)))
            {
                removedGroups.add(phasedRc);
            }
        }

        removedGroups.forEach(x -> mPhasedReadCounters.remove(x));
    }

    public void assignLocalPhaseSets()
    {
        // assign local phase set IDs to all phased variants
        if(!mEnabled)
            return;

        mergeGroups();

        // logPhasedReadCounters();

        for(PhasedReadCounters phasedRc : mPhasedReadCounters.values())
        {
            if(phasedRc.ReadCount + phasedRc.AllocatedReadCount < MIN_READ_COUNT)
                continue;

            int nextLps = mPhaseSetCounter.getNext();
            phasedRc.PositiveReadCounters.forEach(x -> x.addLocalPhaseSet(nextLps, phasedRc.ReadCount, phasedRc.AllocatedReadCount));
        }

        mPhasedReadCounters.clear();
    }

    private static String formRcId(final List<ReadContextCounter> posCounters, final List<ReadContextCounter> negCounters)
    {
        Collections.sort(posCounters, new ReadContextCounterComparator());
        StringJoiner posSj = new StringJoiner("_");
        posCounters.forEach(x -> posSj.add(String.valueOf(x.id())));

        Collections.sort(negCounters, new ReadContextCounterComparator());
        StringJoiner negSj = new StringJoiner("_");
        negCounters.forEach(x -> negSj.add(String.valueOf(x.id())));

        return posSj.toString() + "-" + negSj.toString();
    }

    public static class ReadContextCounterComparator implements Comparator<ReadContextCounter>
    {
        public int compare(final ReadContextCounter first, final ReadContextCounter second)
        {
            return first.id() - second.id();
        }
    }

    private class PhasedReadCounters
    {
        public final int Id;
        public final String RcId;
        public final List<ReadContextCounter> PositiveReadCounters; // supported by the reads
        public final List<ReadContextCounter> NegativeReadCounters; // not supported by the reads

        public int ReadCount; // from uniquely supporting reads
        public double AllocatedReadCount; // allocated from subset groups

        public PhasedReadCounters(
                final int id, final String rcId,
                final List<ReadContextCounter> posCounters, final List<ReadContextCounter> negCounters)
        {
            Id = id;
            RcId = rcId;
            PositiveReadCounters = posCounters;
            NegativeReadCounters = negCounters;
            ReadCount = 1;
            AllocatedReadCount = 0;
        }

        private boolean isExactMatch(final List<ReadContextCounter> posCounters, final List<ReadContextCounter> negCounters)
        {
            if(PositiveReadCounters.size() != posCounters.size() || NegativeReadCounters.size() != negCounters.size())
                return false;

            if(PositiveReadCounters.stream().anyMatch(x -> !posCounters.contains(x)))
                return false;

            if(NegativeReadCounters.stream().anyMatch(x -> !negCounters.contains(x)))
                return false;

            return true;
        }

        private boolean isSubsetOf(final PhasedReadCounters other)
        {
            // returns true if this group is a subset of 'other' is a su
            if(other.PositiveReadCounters.size() < PositiveReadCounters.size())
                return false;

            if(!PositiveReadCounters.stream().allMatch(x -> other.PositiveReadCounters.contains(x)))
                return false;

            // cannot have contradictory negatives
            if(hasAnyOverlap(other.PositiveReadCounters, NegativeReadCounters) || hasAnyOverlap(PositiveReadCounters, other.NegativeReadCounters))
                return false;

            return true;
        }

        private boolean hasAnyOverlap(final List<ReadContextCounter> counters1, final List<ReadContextCounter> counters2)
        {
            return counters1.stream().anyMatch(x -> counters2.contains(x)) || counters2.stream().anyMatch(x -> counters1.contains(x));
        }

        public String toString()
        {
            return String.format("%d: pos(%d) neg(%d) rc(%d) alloc(%.1f)",
                Id, PositiveReadCounters.size(), NegativeReadCounters.size(), ReadCount, AllocatedReadCount);
        }
    }

    private void logPhasedReadCounters()
    {
        Map<ReadContextCounter,Integer> rcIds = Maps.newHashMap();

        for(PhasedReadCounters phasedReadCounters : mPhasedReadCounters.values())
        {
            StringJoiner posVars = new StringJoiner(";");
            StringJoiner posIds = new StringJoiner(";");
            StringJoiner negVars = new StringJoiner(";");
            StringJoiner negIds = new StringJoiner(";");

            for(ReadContextCounter rc : phasedReadCounters.PositiveReadCounters)
            {
                posVars.add(getReadCounterVar(rc));
                posIds.add(String.valueOf(rc.id()));
            }

            for(ReadContextCounter rc : phasedReadCounters.NegativeReadCounters)
            {
                negVars.add(getReadCounterVar(rc));
                negIds.add(String.valueOf(rc.id()));
            }

            SG_LOGGER.debug(String.format("LPS_DATA: %s,%s,%d,%.1f,%s,%s",
                    posIds, negIds, phasedReadCounters.ReadCount, phasedReadCounters.AllocatedReadCount, posVars, negVars));
        }
    }

    private static String getReadCounterVar(final ReadContextCounter rc)
    {
        return String.format("%d_%s>%s", rc.position(), rc.ref(), rc.alt());
    }
}
