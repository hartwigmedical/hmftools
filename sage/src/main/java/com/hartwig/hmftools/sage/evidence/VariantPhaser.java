package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;

public class VariantPhaser
{
    private final PhaseSetCounter mPhaseSetCounter;
    private List<ReadContextCounter> mReadCounters; // only to get consistent indexes

    private boolean mEnabled;

    private final List<PhasedReadCounters> mPhasedReadCounters;

    public static final int MIN_READ_COUNT = 1;

    public VariantPhaser(final PhaseSetCounter phaseSetCounter)
    {
        mPhaseSetCounter = phaseSetCounter;
        mEnabled = false;

        mPhasedReadCounters = Lists.newArrayList();
    }

    public void setEnabled(boolean toggle) { mEnabled = toggle; }
    public boolean enabled() { return mEnabled; }

    public void reset(List<ReadContextCounter> readCounters)
    {
        mPhasedReadCounters.clear();
        mReadCounters = readCounters;
    }

    public void registeredPhasedVariants(final Set<ReadContextCounter> posCounters, final Set<ReadContextCounter> negCounters)
    {
        if(!mEnabled)
            return;

        if(posCounters.isEmpty() || posCounters.size() + negCounters.size() < 2)
            return;

        for(PhasedReadCounters phasedReadCounters : mPhasedReadCounters)
        {
            if(phasedReadCounters.processPhasedCounters(posCounters, negCounters))
                return;
        }

        int id = mPhasedReadCounters.size();
        mPhasedReadCounters.add(new PhasedReadCounters(id, posCounters, negCounters));
    }

    /*
    For each additional row update the set of  candidate LPS using the following rules:  (this will make a maximally collapsed set of candidates)
        If it is a subset of an existing LPS then do nothing
        If it agrees with but extends that LPS then extend the LPS
        If it does not match any overlapping LPS then make a new LPS
    For each LPS, count uniquely matching reads and assigned matching reads (pro rata)
    Filter for LPS support > X
    Remove any uninformative LPS (ie LPS has 1+ variant and is the only LPS that includes that variant).

    Routine:
    - merge the smaller / subsets into the larger groups
    - where multiple options exist, split the counts pro-rata
    - continue until no more merging can be done


    */

    private void mergeGroups()
    {
        List<PhasedReadCounters> removedGroups = Lists.newArrayList();

        for(PhasedReadCounters phasedRc : mPhasedReadCounters)
        {
            if(removedGroups.contains(phasedRc))
                continue;

            // check for super-set groups
            List<PhasedReadCounters> superGroups = mPhasedReadCounters.stream()
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

        // remove any small groups or those with a single variant not present as a +ve in another group
        for(PhasedReadCounters phasedRc : mPhasedReadCounters)
        {
            if(phasedRc.PositiveReadCounters.size() != 1)
                continue;

            if(mPhasedReadCounters.stream()
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

        for(PhasedReadCounters phasedRc : mPhasedReadCounters)
        {
            if(phasedRc.ReadCount + phasedRc.AllocatedReadCount < MIN_READ_COUNT)
                continue;

            int nextLps = mPhaseSetCounter.getNext();
            phasedRc.PositiveReadCounters.forEach(x -> x.addLocalPhaseSet(nextLps, phasedRc.ReadCount, phasedRc.AllocatedReadCount));
        }

        mPhasedReadCounters.clear();
    }

    private class PhasedReadCounters implements Comparable<PhasedReadCounters>
    {
        public final int Id;
        public final Set<ReadContextCounter> PositiveReadCounters; // supported by the reads
        public final Set<ReadContextCounter> NegativeReadCounters; // not supported by the reads

        public int ReadCount; // from uniquely supporting reads
        public double AllocatedReadCount; // allocated from subset groups

        public PhasedReadCounters(final int id, final Set<ReadContextCounter> posCounters, final Set<ReadContextCounter> negCounters)
        {
            Id = id;
            PositiveReadCounters = posCounters;
            NegativeReadCounters = negCounters;
            ReadCount = 1;
            AllocatedReadCount = 0;
        }

        public boolean processPhasedCounters(final Set<ReadContextCounter> posCounters, final Set<ReadContextCounter> negCounters)
        {
            // merge the groups if either +ve or -ve is a subset of the other, with no conflicts
            if(!isExactMatch(posCounters, negCounters))
                return false;

            ReadCount++;
            return true;
        }

        private boolean isExactMatch(final Set<ReadContextCounter> posCounters, final Set<ReadContextCounter> negCounters)
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

        private boolean hasAnyOverlap(final Set<ReadContextCounter> counters1, final Set<ReadContextCounter> counters2)
        {
            return counters1.stream().anyMatch(x -> counters2.contains(x)) || counters2.stream().anyMatch(x -> counters1.contains(x));
        }

        @Override
        public int compareTo(final PhasedReadCounters other)
        {
            return other.ReadCount - ReadCount;
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

        for(PhasedReadCounters phasedReadCounters : mPhasedReadCounters)
        {
            StringJoiner posVars = new StringJoiner(";");
            StringJoiner posIds = new StringJoiner(";");
            StringJoiner negVars = new StringJoiner(";");
            StringJoiner negIds = new StringJoiner(";");

            for(ReadContextCounter rc : phasedReadCounters.PositiveReadCounters)
            {
                posVars.add(getReadCounterVar(rc));
                posIds.add(String.valueOf(getReadCounterId(rcIds, rc)));
            }

            for(ReadContextCounter rc : phasedReadCounters.NegativeReadCounters)
            {
                negVars.add(getReadCounterVar(rc));
                negIds.add(String.valueOf(getReadCounterId(rcIds, rc)));
            }

            SG_LOGGER.debug(String.format("LPS_DATA: %s,%s,%d,%.1f,%s,%s",
                    posIds, negIds, phasedReadCounters.ReadCount, phasedReadCounters.AllocatedReadCount, posVars, negVars));
        }
    }

    private static String getReadCounterVar(final ReadContextCounter rc)
    {
        return String.format("%d_%s>%s", rc.position(), rc.ref(), rc.alt());
    }

    private int getReadCounterId(final Map<ReadContextCounter,Integer> rcIds, final ReadContextCounter rc)
    {
        Integer id = rcIds.get(rc);
        if(id != null)
            return id;

        for(int i = 0; i < mReadCounters.size(); ++i)
        {
            if(rc == mReadCounters.get(i))
            {
                rcIds.put(rc, i);
                return i;
            }
        }

        return -1;
    }
}
