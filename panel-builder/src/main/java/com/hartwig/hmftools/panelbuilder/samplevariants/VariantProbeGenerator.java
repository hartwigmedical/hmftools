package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.Math.min;
import static java.lang.String.format;
import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.panelbuilder.PanelCoverage;
import com.hartwig.hmftools.panelbuilder.Probe;
import com.hartwig.hmftools.panelbuilder.ProbeFactory;
import com.hartwig.hmftools.panelbuilder.ProbeGenerationResult;
import com.hartwig.hmftools.panelbuilder.TargetMetadata;
import com.hartwig.hmftools.panelbuilder.TargetRegion;

import org.jetbrains.annotations.Nullable;

public class VariantProbeGenerator
{
    private record ProbeData(
            String sequence,
            @Nullable ChrBaseRegion start,
            String insert,
            @Nullable ChrBaseRegion end
    )
    {
        List<ChrBaseRegion> regions()
        {
            List<ChrBaseRegion> result = new ArrayList<>();
            if(start != null)
            {
                result.add(start);
            }
            if(end != null)
            {
                result.add(end);
            }
            return result;
        }
    }

    public static ProbeGenerationResult generateMutationProbe(final String chromosome, final int position, final String ref,
            final String alt, final TargetMetadata metadata, final RefGenomeInterface refGenome, final ProbeFactory probeFactory,
            final PanelCoverage coverage)
    {
        ProbeData probeData = generateMutationProbe(refGenome, chromosome, position, ref, alt);
        return createProbe(probeData, metadata, probeFactory, coverage);
    }

    private static ProbeData generateMutationProbe(final RefGenomeInterface refGenome, final String chromosome, final int position,
            final String ref, final String alt)
    {
        int probeLength = PROBE_LENGTH;
        int altLength = alt.length();
        int refLength = ref.length();
        int startLength = probeLength / 2 - altLength / 2;
        int startPos = position - startLength;

        ChrBaseRegion startRegion = new ChrBaseRegion(chromosome, startPos, position - 1);
        String startBases = refGenome.getBaseString(startRegion.chromosome(), startRegion.start(), startRegion.end());
        int endBaseLength = probeLength - startBases.length() - altLength;

        int postPosition = position + refLength;
        ChrBaseRegion endRegion = new ChrBaseRegion(chromosome, postPosition, postPosition + endBaseLength - 1);
        String endBases = refGenome.getBaseString(endRegion.chromosome(), endRegion.start(), endRegion.end());

        String sequence = startBases + alt + endBases;

        if(sequence.length() != probeLength)
        {
            throw new IllegalArgumentException(format("variant(%s:%d %s->%s) invalid sequenceLength(%d): %s",
                    chromosome, position, ref, alt, sequence.length(), sequence));
        }

        return new ProbeData(sequence, startRegion, alt, endRegion);
    }

    public static ProbeGenerationResult generateSvProbe(final String chrStart, final int positionStart, final byte orientStart,
            final String chrEnd, final int positionEnd, final byte orientEnd, final String insertSequence, final TargetMetadata metadata,
            final RefGenomeInterface refGenome, final ProbeFactory probeFactory, final PanelCoverage coverage)
    {
        ProbeData probeData =
                generateSvProbe(refGenome, chrStart, positionStart, orientStart, chrEnd, positionEnd, orientEnd, insertSequence);
        return createProbe(probeData, metadata, probeFactory, coverage);
    }

    private static ProbeData generateSvProbe(final RefGenomeInterface refGenome, final String chrStart, final int positionStart,
            final byte orientStart, final String chrEnd, final int positionEnd, final byte orientEnd, final String insertSequence)
    {
        int probeLength = PROBE_LENGTH;
        int halfProbeLength = probeLength / 2;
        int insSeqLength = insertSequence.length();
        int halfInsSeqLength = insSeqLength / 2;
        int halfNonInsSeqLength = halfProbeLength - halfInsSeqLength;

        // +1 to -1 - do nothing
        // -1 to +1 - a DUP

        // +1 to +1 - start normal, add insert and then reverse compliment of the other side
        // -1 to -1 - alt insert sequence then ref

        ChrBaseRegion startRegion;
        String basesStart;
        ChrBaseRegion endRegion;
        String basesEnd;

        if(orientStart == ORIENT_FWD)
        {
            int probeStart = positionStart - halfNonInsSeqLength + 1;
            startRegion = new ChrBaseRegion(chrStart, probeStart, positionStart);
            basesStart = refGenome.getBaseString(startRegion.chromosome(), startRegion.start(), startRegion.end());

            int endBaseLength = probeLength - basesStart.length() - insSeqLength;

            if(orientEnd == ORIENT_REV)
            {
                endRegion = new ChrBaseRegion(chrEnd, positionEnd, positionEnd + endBaseLength - 1);
                basesEnd = refGenome.getBaseString(endRegion.chromosome(), endRegion.start(), endRegion.end());
            }
            else
            {
                endRegion = new ChrBaseRegion(chrEnd, positionEnd - endBaseLength + 1, positionEnd);
                basesEnd = refGenome.getBaseString(endRegion.chromosome(), endRegion.start(), endRegion.end());
                basesEnd = Nucleotides.reverseComplementBases(basesEnd);
            }
        }
        else
        {
            if(orientEnd == ORIENT_FWD)
            {
                // swap ends and treat as +1/-1
                int probeStart = positionEnd - halfNonInsSeqLength + 1;
                startRegion = new ChrBaseRegion(chrEnd, probeStart, positionEnd);
                basesStart = refGenome.getBaseString(startRegion.chromosome(), startRegion.start(), startRegion.end());

                int endBaseLength = probeLength - basesStart.length() - insSeqLength;
                endRegion = new ChrBaseRegion(chrStart, positionStart, positionStart + endBaseLength - 1);
            }
            else
            {
                // -1/-1 - start with the reversed bases from the end breakend
                startRegion = new ChrBaseRegion(chrEnd, positionEnd, positionEnd + halfNonInsSeqLength - 1);
                basesStart = refGenome.getBaseString(startRegion.chromosome(), startRegion.start(), startRegion.end());
                basesStart = Nucleotides.reverseComplementBases(basesStart);

                int endBaseLength = probeLength - basesStart.length() - insSeqLength;
                endRegion = new ChrBaseRegion(chrStart, positionStart, positionStart + endBaseLength - 1);
            }
            basesEnd = refGenome.getBaseString(endRegion.chromosome(), endRegion.start(), endRegion.end());
        }

        String sequence = basesStart + insertSequence + basesEnd;

        if(sequence.length() != probeLength)
        {
            throw new IllegalArgumentException("Invalid probe");
        }

        return new ProbeData(sequence, startRegion, insertSequence, endRegion);
    }

    public static ProbeGenerationResult generateSglProbe(final String chromosome, final int position, final byte orientation,
            final String insertSequence, final TargetMetadata metadata, final RefGenomeInterface refGenome, final ProbeFactory probeFactory,
            final PanelCoverage coverage)
    {
        ProbeData probeData = generateSglProbe(refGenome, chromosome, position, orientation, insertSequence);
        return createProbe(probeData, metadata, probeFactory, coverage);
    }

    private static ProbeData generateSglProbe(final RefGenomeInterface refGenome, final String chromosome, final int position,
            final byte orientation, final String insertSequence)
    {
        int probeLength = PROBE_LENGTH;
        int halfProbeLength = probeLength / 2;
        int insSeqLength = min(insertSequence.length(), halfProbeLength);
        int refBaseLength = probeLength - insSeqLength;

        // +1 take as-is
        // -1 to +1 - a DUP

        // +1 to +1 - start normal, add insert and then reverse compliment of the other side
        // -1 to -1 - alt insert sequence then ref

        String sequence;

        ChrBaseRegion startRegion = null;
        ChrBaseRegion endRegion = null;
        String insert;

        if(orientation == ORIENT_FWD)
        {
            int probeStart = position - refBaseLength + 1;
            startRegion = new ChrBaseRegion(chromosome, probeStart, position);
            String refBases = refGenome.getBaseString(startRegion.chromosome(), startRegion.start(), startRegion.end());
            insert = insertSequence.substring(0, insSeqLength);
            sequence = refBases + insert;
        }
        else
        {
            insert = insertSequence.substring(insertSequence.length() - insSeqLength);
            endRegion = new ChrBaseRegion(chromosome, position, position + refBaseLength - 1);
            String refBases = refGenome.getBaseString(endRegion.chromosome(), endRegion.start(), endRegion.end());
            sequence = insert + refBases;
        }

        if(sequence.length() != probeLength)
        {
            throw new IllegalArgumentException("Invalid probe");
        }

        return new ProbeData(sequence, startRegion, insert, endRegion);
    }

    private static ProbeGenerationResult createProbe(final ProbeData probeData, final TargetMetadata metadata,
            final ProbeFactory probeFactory, final PanelCoverage coverage)
    {
        List<TargetRegion> targetRegions = probeData.regions().stream()
                .map(region -> new TargetRegion(region, metadata))
                .toList();

        // If there's an alt sequence then always produce a probe to cover it.
        boolean covered = probeData.insert().isEmpty() && probeData.regions().stream().allMatch(coverage::isCovered);
        if(covered)
        {
            return new ProbeGenerationResult(emptyList(), targetRegions, emptyList(), emptyList());
        }
        else
        {
            Probe probe = probeFactory.createProbeFromSequence(probeData.sequence(), metadata).orElseThrow();
            return new ProbeGenerationResult(List.of(probe), targetRegions, targetRegions, emptyList());
        }
    }
}
