package com.hartwig.hmftools.tars.liftback.rescue;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

// Three views of an annotated-junction set: exact membership, by-intron-start, and by-intron-end.
// All three are built from the same Set<ChrIntron> in one pass.
public class AnnotatedJunctionIndex
{
    private final Set<ChrIntron> mJunctions;
    private final Map<ChrPos, List<ChrIntron>> mByStart;
    private final Map<ChrPos, List<ChrIntron>> mByEnd;

    public AnnotatedJunctionIndex(final Set<ChrIntron> junctions)
    {
        mJunctions = junctions != null ? junctions : new HashSet<>();
        mByStart = new HashMap<>();
        mByEnd = new HashMap<>();
        for(ChrIntron intron : mJunctions)
        {
            mByStart.computeIfAbsent(new ChrPos(intron.Chromosome, intron.IntronStart),
                    k -> new ArrayList<>()).add(intron);
            mByEnd.computeIfAbsent(new ChrPos(intron.Chromosome, intron.IntronEnd),
                    k -> new ArrayList<>()).add(intron);
        }
    }

    public boolean contains(final ChrIntron intron)
    {
        return mJunctions.contains(intron);
    }

    public int size()
    {
        return mJunctions.size();
    }

    public List<ChrIntron> introByStart(final String chromosome, final int intronStart)
    {
        return mByStart.getOrDefault(new ChrPos(chromosome, intronStart), Collections.emptyList());
    }

    public List<ChrIntron> introByEnd(final String chromosome, final int intronEnd)
    {
        return mByEnd.getOrDefault(new ChrPos(chromosome, intronEnd), Collections.emptyList());
    }

    private static final class ChrPos
    {
        final String Chromosome;
        final int Position;

        ChrPos(final String chromosome, final int position)
        {
            Chromosome = chromosome;
            Position = position;
        }

        @Override
        public boolean equals(final Object o)
        {
            if(this == o) return true;
            if(!(o instanceof ChrPos)) return false;
            final ChrPos other = (ChrPos) o;
            return Position == other.Position && Chromosome.equals(other.Chromosome);
        }

        @Override
        public int hashCode()
        {
            return Chromosome.hashCode() * 31 + Position;
        }
    }
}
