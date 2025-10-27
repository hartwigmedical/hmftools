package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.genome.gc.GcCalcs.calcGcPercent;
import static com.hartwig.hmftools.panelbuilder.Utils.isDnaSequenceNormal;

import java.util.Optional;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Handles creation of individual probes.
public class ProbeFactory
{
    private final RefGenomeInterface mRefGenome;

    public ProbeFactory(final RefGenomeInterface refGenome)
    {
        mRefGenome = refGenome;
    }

    // Creates a probe from its region(s)/sequence.
    // Returns empty optional if it's not valid to create a probe at that location.
    public Optional<Probe> createProbe(final SequenceDefinition definition, final TargetMetadata metadata)
    {
        // Only check properties which are inconvenient for the caller to check in advance.
        // Everything else is expected to be checked by the caller and will generate an exception.
        boolean regionsValid = definition.regions().stream().allMatch(this::isRegionValid);
        if(!regionsValid)
        {
            // TODO: should probably throw here - why are we trying to create a probe from a bad region?
            return Optional.empty();
        }

        String sequence = buildSequence(definition);
        boolean sequenceValid = sequence.length() == definition.baseLength() && isDnaSequenceNormal(sequence);
        if(!sequenceValid)
        {
            // TODO? maybe throw here. surely very rare to want to create a probe in an unsequenceable region
            return Optional.empty();
        }

        double gcContent = calcGcPercent(sequence);

        // Quality score is calculated later, in batches, because alignment is much faster when batched.

        return Optional.of(new Probe(definition, sequence, metadata, null, null, null, gcContent));
    }

    private boolean isRegionValid(final ChrBaseRegion region)
    {
        return region.hasValidPositions() && region.end() <= mRefGenome.getChromosomeLength(region.chromosome());
    }

    private String buildSequence(final SequenceDefinition definition)
    {
        String start = definition.startRegion() == null ? "" : getSequence(definition.startRegion());
        if(definition.startOrientation() == Orientation.REVERSE)
        {
            start = reverseComplementBases(start);
        }
        String insert = definition.insertSequence() == null ? "" : definition.insertSequence();
        String end = definition.endRegion() == null ? "" : getSequence(definition.endRegion());
        if(definition.endOrientation() == Orientation.REVERSE)
        {
            end = reverseComplementBases(end);
        }
        return start + insert + end;
    }

    private String getSequence(final ChrBaseRegion region)
    {
        return mRefGenome.getBaseString(region.chromosome(), region.start(), region.end());
    }
}
