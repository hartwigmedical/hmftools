package com.hartwig.hmftools.isofox;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class IsofoxConstants
{
    public static final int DEFAULT_MAX_FRAGMENT_SIZE = 550;

    // min number of fragments to sample when calculating fragment length distribution for calculation of expected fragment counts
    public static final int DEFAULT_FRAG_LENGTH_MIN_COUNT = 1000000;

    public static final double GC_RATIO_BUCKET = 0.01;

    public static final short DEFAULT_SINGLE_MAP_QUALITY = 255;
    public static short SINGLE_MAP_QUALITY = DEFAULT_SINGLE_MAP_QUALITY;
    public static short MULTI_MAP_QUALITY_THRESHOLD = 3; // multi-mapped fragments are given map quals of 3 or lower

    public static final int ENRICHED_GENE_BUFFER = 100000;

    public static final int MAX_NOVEL_SJ_DISTANCE = 500000; // beyond which a fragment will be considered chimeric

    public static final double MAX_GENE_PERC_CONTRIBUTION = 0.01;

    // LINC00486
    public static final ChrBaseRegion EXCLUDED_REGION_1_REF_37 = new ChrBaseRegion("2", 33141260, 33141700);
    public static final ChrBaseRegion EXCLUDED_REGION_1_REF_19 = new ChrBaseRegion("chr2", 33141260, 33141700);
    public static final ChrBaseRegion EXCLUDED_REGION_1_REF_38 = new ChrBaseRegion("chr2", 32916190, 32916630);

    public static void populateEnrichedGeneIds(final List<String> geneIds, final RefGenomeVersion version)
    {
        if(version == RefGenomeVersion.V38)
        {
            geneIds.add("ENSG00000276168");
            geneIds.add("ENSG00000274012");
            geneIds.add("ENSG00000278771");
            geneIds.add("ENSG00000263740");
            geneIds.add("ENSG00000283293");
            geneIds.add("ENSG00000265735");
        }
        else
        {
            geneIds.add("ENSG00000265150");
            geneIds.add("ENSG00000258486");
            geneIds.add("ENSG00000202198");
            geneIds.add("ENSG00000266037");
            geneIds.add("ENSG00000263740");
            geneIds.add("ENSG00000265735");
        }
    }

    public static final List<String> ENRICHED_GENE_CHROMOSOMES = Lists.newArrayList("14", "3", "6", "9");
}
