package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;

public class VariantPhaser
{
    private final PhaseSetCounter mPhaseSetCounter;

    private boolean mEnabled;

    private final List<PhasedReadCounters> mPhasedReadCounters;

    public VariantPhaser(final PhaseSetCounter phaseSetCounter)
    {
        mPhaseSetCounter = phaseSetCounter;
        mEnabled = false;

        mPhasedReadCounters = Lists.newArrayList();
    }

    public void setEnabled(boolean toggle) { mEnabled = toggle; }
    public boolean enabled() { return mEnabled; }

    public void reset()
    {
        mPhasedReadCounters.clear();
    }

    public void registeredPhasedVariants(final Set<ReadContextCounter> posCounters, final Set<ReadContextCounter> negCounters)
    {
        if(!mEnabled)
            return;

        if(posCounters.isEmpty() || posCounters.size() + negCounters.size() < 2)
            return;

        for(PhasedReadCounters phasedReadCounters : mPhasedReadCounters)
        {
            if(phasedReadCounters.processPhasedCounters(posCounters, negCounters, 1))
                return;
        }

        mPhasedReadCounters.add(new PhasedReadCounters(posCounters, negCounters));
    }

    public void assignLocalPhaseSets()
    {
        // assign local phase set IDs to all phased variants
        if(!mEnabled)
            return;

        // re-test merge conditions
        for(int i = 0; i < mPhasedReadCounters.size() - 1; ++i)
        {
            PhasedReadCounters phasedRc1 = mPhasedReadCounters.get(i);

            int j = i + 1;
            while(j < mPhasedReadCounters.size())
            {
                PhasedReadCounters phasedRc2 = mPhasedReadCounters.get(j);

                if(phasedRc1.testMerge(phasedRc2))
                    mPhasedReadCounters.remove(j);
                else
                    ++j;
            }
        }

        // logPhasedReadCounters();

        for(PhasedReadCounters phasedReadCounters : mPhasedReadCounters)
        {
            int nextLps = mPhaseSetCounter.getNext();

            // phasedReadCounters.ReadCounters.forEach(x -> x.setLocalPhaseSet(nextLps));
            // phasedReadCounters.ReadCounters.forEach(x -> x.setLpsReadCount(phasedReadCounters.ReadCount));
        }

        mPhasedReadCounters.clear();
    }

    private class PhasedReadCounters
    {
        public final Set<ReadContextCounter> PositiveReadCounters; // supported by the reads
        public final Set<ReadContextCounter> NegativeReadCounters; // not supported by the reads
        public int ReadCount;

        public PhasedReadCounters(final Set<ReadContextCounter> posCounters, final Set<ReadContextCounter> negCounters)
        {
            PositiveReadCounters = posCounters;
            NegativeReadCounters = negCounters;
            ++ReadCount;
        }

        public boolean testMerge(final PhasedReadCounters other)
        {
            return processPhasedCounters(other.PositiveReadCounters, other.NegativeReadCounters, other.ReadCount);
        }

        public boolean processPhasedCounters(
                final Set<ReadContextCounter> posCounters, final Set<ReadContextCounter> negCounters, int readCount)
        {
            // merge the groups if either +ve or -ve is a subset of the other, with no conflicts
            if(!isValidMatch(posCounters, negCounters))
                return false;

            ReadCount += readCount;

            // expand to include the super set
            if(posCounters.size() > PositiveReadCounters.size())
                posCounters.forEach(x -> PositiveReadCounters.add(x));

            if(negCounters.size() > NegativeReadCounters.size())
                negCounters.forEach(x -> NegativeReadCounters.add(x));

            return true;
        }

        private boolean isValidMatch(final Set<ReadContextCounter> posCounters, final Set<ReadContextCounter> negCounters)
        {
            if(!matchesOrSubset(PositiveReadCounters, posCounters))
                return false;

            if(!matchesOrSubset(NegativeReadCounters, negCounters))
                return false;

            // check no negatives exist in the positives and vice versa
            if(hasAnyOverlap(posCounters, NegativeReadCounters) || hasAnyOverlap(PositiveReadCounters, negCounters))
                return false;

            return true;
        }

        private boolean matchesOrSubset(final Set<ReadContextCounter> counters1, final Set<ReadContextCounter> counters2)
        {
            return counters1.stream().allMatch(x -> counters2.contains(x)) || counters2.stream().allMatch(x -> counters1.contains(x));
        }

        private boolean hasAnyOverlap(final Set<ReadContextCounter> counters1, final Set<ReadContextCounter> counters2)
        {
            return counters1.stream().anyMatch(x -> counters2.contains(x)) || counters2.stream().anyMatch(x -> counters1.contains(x));
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

            SG_LOGGER.debug("LPS_DATA: {},{},{},{},{}", posIds, negIds, posVars, negVars, phasedReadCounters.ReadCount);
        }
    }

    private static String getReadCounterVar(final ReadContextCounter rc)
    {
        return String.format("%d_%s>%s", rc.position(), rc.ref(), rc.alt());
    }

    private static int getReadCounterId(final Map<ReadContextCounter,Integer> rcIds, final ReadContextCounter rc)
    {
        Integer id = rcIds.get(rc);
        if(id != null)
            return id;

        rcIds.put(rc, rcIds.size() + 1);
        return rcIds.size();
    }
}
