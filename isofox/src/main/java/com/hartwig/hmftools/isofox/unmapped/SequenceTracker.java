package com.hartwig.hmftools.isofox.unmapped;

import java.util.List;

import com.google.common.collect.Lists;

public class SequenceTracker
{
    private final List<ConsensusSequence> mSequences;
    private final List<UnmappedRead> mAllReads; // ordered by longest SC length descending

    private static final int MAX_SEQUENCE_MERGE_BASE_DIFFS = 2;

    public SequenceTracker()
    {
        mSequences = Lists.newArrayList();
        mAllReads = Lists.newArrayList();
    }

    public void clear()
    {
        mSequences.clear();
        mAllReads.clear();
    }

    public List<ConsensusSequence> getSequences() { return mSequences; }

    public void processRead(final UnmappedRead umRead)
    {
        /*
        Some small notes on consensus algo:
        - Traverse starting from 1st intronic base, voting on each base in turn
        - If multiple high quality bases deviate from sequence then add a 2nd consensus sequence and divide reads amongst those groups.
        - If a single read deviates for more than 2 high qual bases within say 5 bases (?) then add a separate consensus also
        - Check consensus sequences at end for similarity and try to collapse (eg. due to DEL or INS).
        */

        boolean matched = false;

        for(ConsensusSequence sequence : mSequences)
        {
            if(sequence.checkAddReadBases(umRead))
            {
                matched = true;
            }
        }

        if(!matched)
        {
            mSequences.add(new ConsensusSequence(umRead));
        }

        int index = 0;
        while(index < mAllReads.size())
        {
            if(umRead.ScBases.length() > mAllReads.get(index).ScBases.length())
                break;

            ++index;
        }

        mAllReads.add(index, umRead);
    }

    public void reconcileSequences()
    {
        List<UnmappedRead> unassignedReads = Lists.newArrayList(mAllReads);

        for(ConsensusSequence sequence : mSequences)
        {
            sequence.setExactMatchReads();

            sequence.exactMatchReads().forEach(x -> unassignedReads.remove(x));
        }

        for(UnmappedRead umRead : unassignedReads)
        {
            List<ConsensusSequence> candidateSequences = Lists.newArrayList();

            double matchTotal = 0;

            for(ConsensusSequence sequence : mSequences)
            {
                if(sequence.matches(umRead, false))
                {
                    candidateSequences.add(sequence);
                    matchTotal += sequence.assignedReadTotal();
                }
            }

            for(ConsensusSequence sequence : candidateSequences)
            {
                double portion = sequence.assignedReadTotal() / matchTotal;
                sequence.addInexactRead(umRead, portion);
            }
        }

        // merge any sequence with minimal support
        List<ConsensusSequence> mergedSequences = Lists.newArrayList();

        for(int i = 0; i < mSequences.size() - 1; ++i)
        {
            ConsensusSequence seq1 = mSequences.get(i);

            if(mergedSequences.contains(seq1))
                continue;

            for(int j = i + 1; j < mSequences.size(); ++j)
            {
                ConsensusSequence seq2 = mSequences.get(j);

                if(mergedSequences.contains(seq2))
                    continue;

                if(seq1.getBaseDiffCount(seq2.getSequenceString()) > MAX_SEQUENCE_MERGE_BASE_DIFFS)
                    continue;

                int unique1 = uniqueReads(seq1, seq2);
                int unique2 = uniqueReads(seq2, seq1);

                if(unique1 <= 1 || unique2 <= 1)
                {
                    if(unique1 >= unique2)
                    {
                        mergedSequences.add(seq2);
                        seq1.reassignReads(seq2);
                    }
                    else
                    {
                        mergedSequences.add(seq1);
                        seq2.reassignReads(seq1);
                        break;
                    }
                }
            }
        }

        mergedSequences.forEach(x -> mSequences.remove(x));
    }

    private static int uniqueReads(final ConsensusSequence seq, final ConsensusSequence other)
    {
        // count reads in the first sequence which are not in the second
        return (int)seq.allMatchedReads().stream()
                .filter(x -> !other.exactMatchReads().contains(x) && !other.inexactMatchReads().contains(x))
                .count();
    }
}
