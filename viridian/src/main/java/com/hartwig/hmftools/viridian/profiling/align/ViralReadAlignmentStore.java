package com.hartwig.hmftools.viridian.profiling.align;

import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

import com.hartwig.hmftools.viridian.profiling.common.ReadId;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

// A set of read alignments to virus genomes.
// For each read, stores only the best alignment to each virus genome.
// (Rarely, a read can align to the same virus genome multiple times, but these secondary alignments are not
// necessary to keep).
public class ViralReadAlignmentStore
{
    private final Map<VirusGenome, Map<ReadId, ViralReadAlignment>> mByGenome;
    private final Map<ReadId, Map<VirusGenome, ViralReadAlignment>> mByRead;

    public ViralReadAlignmentStore(Collection<ViralReadAlignment> alignments)
    {
        mByGenome = new HashMap<>();
        mByRead = new LinkedHashMap<>();
        for(ViralReadAlignment alignment : alignments)
        {
            ViralReadAlignment replaced = mByGenome
                    .computeIfAbsent(alignment.genome(), k -> new LinkedHashMap<>())
                    .put(alignment.readId(), alignment);
            if(replaced != null)
            {
                throw new IllegalArgumentException(
                        "Repeat alignment of a read on virus genome: " + alignment.genome().name());
            }
            mByRead.computeIfAbsent(alignment.readId(), k -> new HashMap<>()).put(alignment.genome(), alignment);
        }
    }

    // Virus genome -> read -> read's best alignment on virus genome.
    public Map<VirusGenome, Map<ReadId, ViralReadAlignment>> byGenome()
    {
        return mByGenome;
    }

    // Read -> virus genome -> read's best alignment on virus genome.
    public Map<ReadId, Map<VirusGenome, ViralReadAlignment>> byRead()
    {
        return mByRead;
    }
}
