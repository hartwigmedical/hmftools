package com.hartwig.hmftools.pave.transval;

import java.util.HashSet;
import java.util.Set;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import org.jetbrains.annotations.NotNull;

class Duplication extends ProteinVariant
{
    public Duplication(@NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            @NotNull final AminoAcidRange refRange)
    {
        super(gene, transcript, aminoAcidSequence, refRange.startPosition(), refRange.length());
    }

    @Override
    TransvalVariant calculateVariant(final RefGenomeInterface genome)
    {
        ChangeContext changeContext = getChangeContext(genome);
        String duplicated = changeContext.basesForProteinChange(positionOfFirstAlteredCodon(), RefLength, Transcript.posStrand()).segmentThatIsModified();
        String baseAtChangeLocation = "";
        Set<TransvalHotspot> hotspots = new HashSet<>();
//        hotspots.add(changeContext.hotspot())
        return new TransvalVariant(
                Transcript,
                Gene.Chromosome,
                false,
                hotspots);
    }
}
