package com.hartwig.hmftools.tars.liftback.supplementary;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Three views of an annotated-junction set: exact membership, by-intron-start, and by-intron-end.
// All three are built from the same Set<ChrBaseRegion> in one pass.
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
            mByStart.computeIfAbsent(new BasePosition(intron.Chromosome, intron.start()),
                    k -> new ArrayList<>()).add(intron);
            mByEnd.computeIfAbsent(new BasePosition(intron.Chromosome, intron.end()),
                    k -> new ArrayList<>()).add(intron);
        }
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
