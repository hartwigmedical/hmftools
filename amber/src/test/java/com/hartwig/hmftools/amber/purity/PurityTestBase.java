package com.hartwig.hmftools.amber.purity;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.amber.PositionEvidence;
import com.hartwig.hmftools.common.amber.AmberBase;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

public class PurityTestBase
{
    PositionEvidence evidenceWithDepthAndAltCount(HumanChromosome chromosome, int position, int depth, int altCount)
    {
        Preconditions.checkArgument(altCount <= depth);
        final PositionEvidence positionEvidence = new PositionEvidence(V38.versionedChromosome(chromosome), position, "A", "T");
        positionEvidence.ReadDepth = depth;
        positionEvidence.RefSupport = depth - altCount;
        positionEvidence.AltSupport = altCount;
        return positionEvidence;
    }

    PositionEvidence evidenceWithDepthAndAltCount(AmberBase ref, AmberBase alt, int position, int depth, int altCount)
    {
        final PositionEvidence positionEvidence = new PositionEvidence(V38.versionedChromosome(_1), position, ref.name(), alt.name());
        positionEvidence.ReadDepth = depth;
        positionEvidence.RefSupport = depth - altCount;
        positionEvidence.AltSupport = altCount;
        return positionEvidence;
    }

    PositionEvidence evidenceWithDepthAndAltCount(int depth, int altCount)
    {
        return evidenceWithDepthAndAltCount(_1, 1234, depth, altCount);
    }

    PositionEvidence evidenceWithDepth(int depth)
    {
        final PositionEvidence positionEvidence = new PositionEvidence("1", 1234, "A", "T");
        positionEvidence.ReadDepth = depth;
        positionEvidence.RefSupport = depth;
        return positionEvidence;
    }
}
