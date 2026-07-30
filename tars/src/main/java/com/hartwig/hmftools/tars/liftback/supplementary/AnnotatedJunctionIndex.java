package com.hartwig.hmftools.tars.liftback.supplementary;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.common.ContigEntry;

// Three views of the annotated junctions: exact membership, by-intron-start and by-intron-end, all built in one pass.
// The sidecar carries exons, so the junctions are the gaps between each transcript's adjacent exon spans.
public class AnnotatedJunctionIndex
{
    private final Set<ChrBaseRegion> mJunctions;
    private final Map<BasePosition, List<ChrBaseRegion>> mByStart;
    private final Map<BasePosition, List<ChrBaseRegion>> mByEnd;

    public AnnotatedJunctionIndex(final Set<ChrBaseRegion> junctions)
    {
        mJunctions = junctions != null ? junctions : new HashSet<>();
        mByStart = new HashMap<>();
        mByEnd = new HashMap<>();
        for(ChrBaseRegion intron : mJunctions)
        {
            mByStart.computeIfAbsent(
                    new BasePosition(intron.Chromosome, intron.start()),
                    k -> new ArrayList<>()).add(intron);
            mByEnd.computeIfAbsent(
                    new BasePosition(intron.Chromosome, intron.end()),
                    k -> new ArrayList<>()).add(intron);
        }
    }

    // Ensembl reaches liftback as the sidecar's exon spans, so the junctions are the gaps between each transcript's
    // adjacent exons. A single-exon transcript contributes none, so nothing is missed by deriving them this way.
    public static AnnotatedJunctionIndex fromContigEntries(final List<ContigEntry> entries)
    {
        Set<ChrBaseRegion> junctions = new HashSet<>();

        for(ContigEntry entry : entries)
        {
            List<BaseRegion> exons = new ArrayList<>(entry.exonSpans());
            if(exons.size() < 2)
            {
                continue;
            }

            exons.sort(Comparator.comparingInt(BaseRegion::start));

            for(int i = 0; i < exons.size() - 1; ++i)
            {
                int intronStart = exons.get(i).end() + 1;
                int intronEnd = exons.get(i + 1).start() - 1;
                if(intronEnd >= intronStart)   // abutting exons leave no intron
                {
                    junctions.add(new ChrBaseRegion(entry.chromosome(), intronStart, intronEnd));
                }
            }
        }

        return new AnnotatedJunctionIndex(junctions);
    }

    public boolean contains(final ChrBaseRegion intron)
    {
        return mJunctions.contains(intron);
    }

    public int size()
    {
        return mJunctions.size();
    }

    public List<ChrBaseRegion> introByStart(final String chromosome, final int intronStart)
    {
        return mByStart.getOrDefault(new BasePosition(chromosome, intronStart), Collections.emptyList());
    }

    public List<ChrBaseRegion> introByEnd(final String chromosome, final int intronEnd)
    {
        return mByEnd.getOrDefault(new BasePosition(chromosome, intronEnd), Collections.emptyList());
    }
}
