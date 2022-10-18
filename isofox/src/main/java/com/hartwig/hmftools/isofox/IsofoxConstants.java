package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.common.immune.ImmuneRegions.getIgRegions;

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

    public static final int ENRICHED_GENE_BUFFER = 1000;

    public static final int MAX_NOVEL_SJ_DISTANCE = 500000; // beyond which a fragment will be considered chimeric

    public static final double MAX_GENE_PERC_CONTRIBUTION = 0.01;

    public static void populateEnrichedGeneIds(final List<String> geneIds, final RefGenomeVersion version)
    {
        if(version == RefGenomeVersion.V38)
        {
            geneIds.add("ENSG00000276168"); // RN7SL1
            geneIds.add("ENSG00000274012"); // RN7SL2
            geneIds.add("ENSG00000278771"); // RN7SL3
            geneIds.add("ENSG00000263740"); // RN7SL4P
            geneIds.add("ENSG00000283293"); // RN7SK
            geneIds.add("ENSG00000265735"); // RN7SL5P
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

    // process genes with longer expected processing times first
    public static final List<String> PRIORITISED_CHROMOSOMES = Lists.newArrayList("14", "22", "2", "3", "6", "9", "16", "17");

    public static void populateImmuneRegions(final List<ChrBaseRegion> regions, final RefGenomeVersion version)
    {
        regions.addAll(getIgRegions(version));

        // also add in the HLA genes region
        if(version == RefGenomeVersion.V38)
            regions.add(new ChrBaseRegion("chr6", 29722775, 33143325));
        else
            regions.add(new ChrBaseRegion("6", 29690552, 33111102));
    }
}
