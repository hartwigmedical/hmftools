package com.hartwig.hmftools.redux.ms_sites;

import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import static htsjdk.samtools.util.SequenceUtil.N;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.ListIterator;
import java.util.concurrent.Callable;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.redux.jitter.JitterConstants;
import com.hartwig.hmftools.redux.jitter.MicrosatelliteSite;

import org.apache.commons.lang3.mutable.MutableInt;

public class ChromosomeSiteFinder implements Callable<Void>
{
    private final String mChromosome;
    private final int mChromosomeLength;
    private final int mPartitionSize;
    private final RefGenomeInterface mRefGenome;
    final Consumer<MicrosatelliteSite> mConsumer;

    public ChromosomeSiteFinder(
            final RefGenomeInterface refGenome, final String chromosome, final int chromosomeLength, final int partitionSize,
            final Consumer<MicrosatelliteSite> consumer)
    {
        mRefGenome = refGenome;
        mChromosome = chromosome;
        mChromosomeLength = chromosomeLength;
        mPartitionSize = partitionSize;
        mConsumer = consumer;
    }

    @Override
    public Void call()
    {
        findSites();
        return null;
    }

    public void findSites()
    {
        // pending candidate, we do not accept a candidate when it is completed. We want to avoid
        // candidates that are too close to each other. They are dropped if too close.
        List<MicrosatelliteSite> pendingMicrosatellies = new ArrayList<>();

        // current best candidate
        Candidate bestCandidate = null;

        // candidates
        List<Candidate> currentCandidates = new ArrayList<>();

        for(int start = 1; start <= mChromosomeLength; start = start + mPartitionSize)
        {
            int endInclusive = Math.min(start + mPartitionSize, mChromosomeLength) - 1;
            byte[] seq = mRefGenome.getBases(mChromosome, start, endInclusive);

            for(int i = 0; i < seq.length; ++i)
            {
                byte base = seq[i];

                if(base == N)
                {
                    // if we hit an N, we clear everything
                    currentCandidates.clear();
                    bestCandidate = null;
                    continue;
                }

                // check all current candidates
                ListIterator<Candidate> itr = currentCandidates.listIterator();
                while(itr.hasNext())
                {
                    Candidate c = itr.next();

                    // check if pattern is still valid
                    if(c.nextBase() == base)
                    {
                        // pattern continues
                        c.currentEndIndex++;
                        if(bestCandidate == null || bestCandidate.numFullUnits() < c.numFullUnits())
                        {
                            bestCandidate = c;
                        }
                    }
                    else
                    {
                        c.complete = true;
                        itr.remove();
                    }
                }

                // if the best candidate is completed, we want to check if it is a candidate
                // we want to include
                if(bestCandidate != null && bestCandidate.complete)
                {
                    // check if the best candidate that is completed is a good one
                    int unitRepeatCount = bestCandidate.numFullUnits();

                    if(unitRepeatCount >= JitterConstants.MIN_MICROSAT_UNIT_COUNT)
                    {
                        // NOTE: for microsatellites with a pattern longer than one, only accept full repeats
                        // i.e. ATGATGATGATGAT contains 4 full ATG plus AT at the end, even though AT is partial unit, they are excluded

                        int baseLength = unitRepeatCount * bestCandidate.pattern.length;

                        // this is a microsatellite
                        MicrosatelliteSite microsatelliteSite = new MicrosatelliteSite(
                                mChromosome,
                                bestCandidate.startIndex,
                                bestCandidate.startIndex + baseLength - 1, // change to inclusive
                                bestCandidate.pattern);
                        pendingMicrosatellies.add(microsatelliteSite);

                        // check the panding microsatellites, see if any can be accepted
                        checkPendingMicrosatellites(pendingMicrosatellies);
                    }

                    bestCandidate = null;
                }

                // also start new Candidates at this location
                // all these candidates have the current base as the last base of the pattern
                for(int j = Math.max(i - JitterConstants.MAX_MICROSAT_UNIT_LENGTH + 1, 0); j <= i; ++j)
                {
                    byte[] repeatUnit = Arrays.copyOfRange(seq, j, i + 1);

                    // check that this is a valid repeat unit, i.e. it is not multiple smaller unit
                    if(isValidUnit(repeatUnit))
                    {
                        // also check against existing patterns
                        if(currentCandidates.stream().noneMatch(o -> Arrays.equals(o.pattern, repeatUnit)))
                        {
                            Candidate c = new Candidate(repeatUnit, start + j, start + i + 1);
                            currentCandidates.add(c);
                        }
                    }
                }
            }
        }

        if(pendingMicrosatellies.size() == 1)
        {
            mConsumer.accept(pendingMicrosatellies.get(0));
        }

        RD_LOGGER.info("finished chromosome {}", mChromosome);
    }

    private void checkPendingMicrosatellites(final List<MicrosatelliteSite> pendingMicrosatellies)
    {
        // the aim of this code is to remove microsatellites that are too close to each other eg AAAAATAAAAAAA
        int groupStart = 0;

        // check the pending MS sites, see if any can be accepted
        for(int i = 0; i < pendingMicrosatellies.size() - 1; ++i)
        {
            MicrosatelliteSite ms1 = pendingMicrosatellies.get(i);
            MicrosatelliteSite ms2 = pendingMicrosatellies.get(i + 1);

            if((ms2.referenceStart() - ms1.referenceEnd()) > JitterConstants.MIN_ADJACENT_MICROSAT_DISTANCE)
            {
                // previous group finished, if previous group only has 1 ms, we accept it, otherwise
                // remove them all from the pending list
                if(i == groupStart)
                {
                    // only 1 item, accept this
                    mConsumer.accept(ms1);
                    // microsatelliteCounter.increment();
                    RD_LOGGER.trace("microsatellite: {}", ms1);
                }
                else
                {
                    // ms are too close to each other, remove them
                    for(int j = groupStart; j <= i; ++j)
                    {
                        RD_LOGGER.trace("reject microsatellite as too close to neighbour: {}", pendingMicrosatellies.get(j));
                    }
                }

                // update group start
                groupStart = i + 1;
            }
        }

        if(groupStart > 0)
            pendingMicrosatellies.subList(0, groupStart).clear();
    }

    private static boolean isValidUnit(byte[] unit)
    {
        // check if this is a valid repeat unit, eg ATAT is not because it is actually 2 x AT

        // check each subunit length to see if it is the same one
        for(int subunitLength = 1; subunitLength <= unit.length / 2; ++subunitLength)
        {
            boolean allMatch = true;

            if((unit.length % subunitLength) == 0)
            {
                for(int i = 0; i < subunitLength && allMatch; ++i)
                {
                    byte base = unit[i];
                    for(int j = subunitLength + i; j < unit.length; j += subunitLength)
                    {
                        if(unit[j] != base)
                        {
                            allMatch = false;
                            break;
                        }
                    }
                }
            }
            else
            {
                allMatch = false;
            }

            if(allMatch)
            {
                // found a subunit that matches the whole unit
                return false;
            }
        }
        // did not find a subunit
        return true;
    }
}
