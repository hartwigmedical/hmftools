package com.hartwig.hmftools.tars.liftback;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.common.ContigEntry;

// Sidecar-derived annotation used during liftback: merged exon spans plus exact intron/junction membership.
public final class EnsemblAnnotationIndex
{
    private final Map<String, List<BaseRegion>> mExonsByChromosome;
    private final Set<ChrBaseRegion> mJunctions;
    private final Map<ChrBaseRegion, Integer> mJunctionStrands;

    private EnsemblAnnotationIndex(
            final Map<String, List<BaseRegion>> exonsByChromosome, final Set<ChrBaseRegion> junctions,
            final Map<ChrBaseRegion, Integer> junctionStrands)
    {
        mExonsByChromosome = exonsByChromosome;
        mJunctions = junctions != null ? Set.copyOf(junctions) : Collections.emptySet();
        mJunctionStrands = junctionStrands != null ? Map.copyOf(junctionStrands) : Collections.emptyMap();
    }

    public boolean containsExon(final String chromosome, final int pos)
    {
        List<BaseRegion> exons = mExonsByChromosome.get(chromosome);
        if(exons == null)
            return false;

        // merged spans are non-overlapping, so the last exon starting at or before pos is the only candidate.
        BaseRegion exon = exons.get(BaseRegion.binarySearch(pos, exons));
        return exon.containsPosition(pos);
    }

    public boolean containsJunction(final ChrBaseRegion intron)
    {
        return mJunctions.contains(intron);
    }

    public int junctionStrand(final ChrBaseRegion intron)
    {
        return mJunctionStrands.getOrDefault(intron, 0);
    }

    public int junctionCount()
    {
        return mJunctions.size();
    }

    public static EnsemblAnnotationIndex fromJunctions(final Set<ChrBaseRegion> junctions)
    {
        return new EnsemblAnnotationIndex(new HashMap<>(), junctions, Collections.emptyMap());
    }

    public static EnsemblAnnotationIndex fromContigEntries(final List<ContigEntry> entries)
    {
        Map<String, List<BaseRegion>> exonsByChromosome = new HashMap<>();
        Set<ChrBaseRegion> junctions = new HashSet<>();
        Map<ChrBaseRegion, Integer> junctionStrands = new HashMap<>();

        for(ContigEntry entry : entries)
        {
            addEntry(entry, exonsByChromosome, junctions, junctionStrands);
        }

        exonsByChromosome.values().forEach(EnsemblAnnotationIndex::mergeSpans);
        return new EnsemblAnnotationIndex(exonsByChromosome, junctions, junctionStrands);
    }

    private static void addEntry(
            final ContigEntry entry, final Map<String, List<BaseRegion>> exonsByChromosome,
            final Set<ChrBaseRegion> junctions, final Map<ChrBaseRegion, Integer> junctionStrands)
    {
        List<BaseRegion> exons = copyExons(entry.exonSpans());
        exonsByChromosome.computeIfAbsent(entry.chromosome(), k -> new ArrayList<>()).addAll(exons);
        addJunctions(junctions, junctionStrands, entry.chromosome(), entry.strand(), exons);
    }

    private static List<BaseRegion> copyExons(final List<BaseRegion> exons)
    {
        List<BaseRegion> copy = new ArrayList<>(exons.size());
        for(BaseRegion exon : exons)
        {
            copy.add(new BaseRegion(exon.start(), exon.end()));
        }
        return copy;
    }

    private static void addJunctions(
            final Set<ChrBaseRegion> junctions, final Map<ChrBaseRegion, Integer> junctionStrands,
            final String chromosome, final int strand, final List<BaseRegion> exons)
    {
        if(exons.size() < 2)
        {
            return;
        }

        exons.sort(Comparator.comparingInt(BaseRegion::start));
        for(int i = 0; i < exons.size() - 1; ++i)
        {
            int intronStart = exons.get(i).end() + 1;
            int intronEnd = exons.get(i + 1).start() - 1;
            if(intronEnd >= intronStart)
            {
                ChrBaseRegion intron = new ChrBaseRegion(chromosome, intronStart, intronEnd);
                junctions.add(intron);
                junctionStrands.merge(intron, strand, (existing, incoming) -> existing.equals(incoming) ? existing : 0);
            }
        }
    }

    // Union-merge overlapping/adjacent spans into sorted, non-overlapping ranges. Not BaseRegion.checkMergeOverlaps:
    // it sets the merged end from the next span rather than the max, so an exon nested inside an earlier-starting one
    // shrinks the merged span - and exons nest across isoforms.
    private static void mergeSpans(final List<BaseRegion> regions)
    {
        regions.sort(Comparator.comparingInt(BaseRegion::start));

        int index = 0;
        while(index + 1 < regions.size())
        {
            BaseRegion current = regions.get(index);
            BaseRegion next = regions.get(index + 1);
            if(next.start() <= current.end() + 1)
            {
                current.setEnd(Math.max(current.end(), next.end()));
                regions.remove(index + 1);
            }
            else
            {
                ++index;
            }
        }
    }
}
