package com.hartwig.hmftools.lilac.qc;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.DELIM;
import static com.hartwig.hmftools.lilac.LilacUtils.listMax;
import static com.hartwig.hmftools.lilac.LilacUtils.listMin;

import java.util.Comparator;
import java.util.List;

import com.hartwig.hmftools.lilac.seq.SequenceCount;
import org.apache.commons.math3.util.Pair;

public class Haplotype
{
    public final int StartLocus;
    public final int EndLocus;
    public final int PhasedFragmentCount;
    public final String Haplotype;

    private int mMatchingFragmentCount;

    public Haplotype(final int startLocus, final int endLocus, final int phasingFragments, final String haplotype)
    {
        StartLocus = startLocus;
        EndLocus = endLocus;
        PhasedFragmentCount = phasingFragments;
        Haplotype = haplotype;
        mMatchingFragmentCount = 0;
    }

    public int matchingFragmentCount() { return mMatchingFragmentCount; }
    public void addMatchingFragmentCount() { ++mMatchingFragmentCount; }

    public static Haplotype fromString(final String line)
    {
        String[] items = line.split(DELIM, -1);
        int startLocus = Integer.parseInt(items[0]);
        String haplotype = items[1];
        int endLocus = startLocus + haplotype.length() - 1;

        return new Haplotype(startLocus, endLocus, 0, haplotype);
    }

    public static Haplotype create(final List<Integer> aminoAcidIndices, final Pair<String,Integer> evidence, final SequenceCount aminoAcidCount)
    {
        int startLocus = listMin(aminoAcidIndices);
        int endLocus = listMax(aminoAcidIndices);
        String sparseHaplotype = evidence.getFirst();
        String complexHaplotype = "";

        for(int locus = startLocus; locus <= endLocus; ++locus)
        {
            if(aminoAcidIndices.contains(locus))
                complexHaplotype += sparseHaplotype.charAt(aminoAcidIndices.indexOf(locus));
            else
                complexHaplotype += aminoAcidCount.getMaxCountSequence(locus);
        }

        return new Haplotype(startLocus, endLocus, evidence.getSecond(), complexHaplotype);
    }

    public boolean contains(final Haplotype unmatched)
    {
        if (unmatched.StartLocus >= StartLocus && unmatched.EndLocus <= EndLocus)
        {
            int start = max(unmatched.StartLocus, StartLocus);
            int end = min(unmatched.EndLocus, EndLocus);

            for(int locus = start; locus <= end; ++locus)
            {
                if(charAt(locus) != unmatched.charAt(locus))
                    return false;
            }

            return true;
        }

        return false;
    }

    private char charAt(int locus)
    {
        int index = locus - StartLocus;
        return Haplotype.charAt(index);
    }

    public String sequence(int locus)
    {
        if(!positionWithin(locus, StartLocus, EndLocus))
            return "";

        int index = locus - StartLocus;

        if(index < 0 || index >= Haplotype.length())
        {
            LL_LOGGER.error("haplotype({}) invalid sequence request at locus({})",
                    toString(), locus);
            return "";
        }

        return String.valueOf(Haplotype.charAt(index));
    }

    public String toString()
    {
        return String.format("loci(%d - %d) frags(supporting=%d phased=%d) haplotype=%s",
                StartLocus, EndLocus, mMatchingFragmentCount, PhasedFragmentCount, Haplotype);
    }

    public static class HaplotypeFragmentSorter implements Comparator<Haplotype>
    {
        public int compare(final Haplotype first, final Haplotype second)
        {
            // by fragments descending
            int compare = first.matchingFragmentCount() - second.matchingFragmentCount();
            if(compare != 0)
                return compare < 0 ? 1 : -1;

            return 0;
        }
    }

    public static class HaplotypeStartLocusSorter implements Comparator<Haplotype>
    {
        public int compare(final Haplotype first, final Haplotype second)
        {
            int compare = first.StartLocus - second.StartLocus;
            if(compare != 0)
                return compare > 0 ? 1 : -1; // lower start first

            return 0;
        }

    }

}
