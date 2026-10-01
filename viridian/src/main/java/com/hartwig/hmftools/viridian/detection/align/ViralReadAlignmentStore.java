package com.hartwig.hmftools.viridian.detection.align;

import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

import com.hartwig.hmftools.viridian.detection.common.ReadId;
import com.hartwig.hmftools.viridian.reference.ViralContig;

// A set of read alignments to virus genomes.
// For each read, stores on the best alignment to each contig.
// (Rarely, a read can align to the same contig multiple times, but these secondary alignments are not necessary to keep).
public class ViralReadAlignmentStore
{
    private final Map<ViralContig, Map<ReadId, ViralReadAlignment>> mByContig;
    private final Map<ReadId, Map<ViralContig, ViralReadAlignment>> mByRead;

    public ViralReadAlignmentStore(Collection<ViralReadAlignment> alignments)
    {
        mByContig = new HashMap<>();
        mByRead = new LinkedHashMap<>();
        for(ViralReadAlignment alignment : alignments)
        {
            ViralReadAlignment replaced = mByContig
                    .computeIfAbsent(alignment.contig(), k -> new LinkedHashMap<>())
                    .put(alignment.readId(), alignment);
            if(replaced != null)
            {
                throw new IllegalArgumentException("Repeat alignment of a read on contig: " + alignment.contig().name());
            }
            mByRead.computeIfAbsent(alignment.readId(), k -> new HashMap<>()).put(alignment.contig(), alignment);
        }
    }

    // Contig -> read -> read's best alignment on contig.
    public Map<ViralContig, Map<ReadId, ViralReadAlignment>> byContig()
    {
        return mByContig;
    }

    // Read -> contig -> read's best alignment on contig.
    public Map<ReadId, Map<ViralContig, ViralReadAlignment>> byRead()
    {
        return mByRead;
    }

    public int readCount()
    {
        return mByRead.size();
    }
}
