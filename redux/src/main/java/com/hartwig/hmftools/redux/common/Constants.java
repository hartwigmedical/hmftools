package com.hartwig.hmftools.redux.common;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.MT_CHR_V37;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.MT_CHR_V38;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class Constants
{
    public static final String FILE_ID = "redux";

    public static final int DEFAULT_READ_LENGTH = 151;

    // UMIs
    public static final int DEFAULT_MAX_UMI_BASE_DIFF = 1;
    public static final int MAX_IMBALANCED_UMI_BASE_DIFF = 4;
    public static final int MAX_IMBALANCED_UMI_COUNT = 25;
    public static final int MIN_POLYG_UMI_TAIL_LENGTH = 2;
    public static final int MAX_UMI_BASE_DIFF_JITTER_COLLAPSE = 1;

    public static final char DEFAULT_DUPLEX_UMI_DELIM = '_';

    public static final String CONSENSUS_PREFIX = "CNS_";

    // read unmapping
    public static final int UNMAP_MAX_NON_OVERLAPPING_BASES = 10;
    public static final int UNMAP_MIN_SOFT_CLIP = 20;
    public static final int UNMAP_MIN_HIGH_DEPTH = 1000;
    public static final int UNMAP_CHIMERIC_FRAGMENT_LENGTH_MAX = 1000;

    public static final int SUPP_ALIGNMENT_SCORE_MIN = 30;

    // consensus building
    public static int CONSENSUS_MAX_DEPTH = 100;

    public static final List<String> EXPECTED_VALID_CONTIGS = Lists.newArrayListWithExpectedSize(30);

    public static void loadExpectedValidContigs(final RefGenomeVersion refGenomeVersion)
    {
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            EXPECTED_VALID_CONTIGS.add(refGenomeVersion.versionedChromosome(chromosome.toString()));
        }

        EXPECTED_VALID_CONTIGS.add(refGenomeVersion.is37() ? MT_CHR_V37 : MT_CHR_V38);
        EXPECTED_VALID_CONTIGS.add(refGenomeVersion.versionedChromosome("EBV"));
    }

}
