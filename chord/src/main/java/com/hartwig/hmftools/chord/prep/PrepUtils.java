package com.hartwig.hmftools.chord.prep;

import static com.hartwig.hmftools.chord.ChordConstants.CHORD_LOGGER;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import htsjdk.variant.variantcontext.VariantContext;

public class PrepUtils
{
    public static void checkRefGenomeVersion(RefGenomeSource refGenome, VariantContext variantContext)
    {
        boolean hasChrPrefix = variantContext.getContig().toLowerCase().startsWith(HumanChromosome.CHR_PREFIX);

        RefGenomeVersion expectedVersion = RefGenomeSource.deriveRefGenomeVersion(refGenome);
        RefGenomeVersion likelyActualVersion = hasChrPrefix ? RefGenomeVersion.V38 : RefGenomeVersion.V37;

        if(expectedVersion != likelyActualVersion)
        {
            CHORD_LOGGER.warn("Expected ref genome version {} but input file is likely {} (detected chromosome name {} prefix '{}')",
                    expectedVersion,
                    likelyActualVersion,
                    hasChrPrefix ? "with" : "without",
                    HumanChromosome.CHR_PREFIX
            );
        }
    }
}
