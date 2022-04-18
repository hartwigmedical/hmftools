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

    public static final List<String> ENRICHED_GENE_CHROMOSOMES = Lists.newArrayList("14", "22", "2", "3", "6", "9");

    // HLA region: chr 6: 29690552-33111102
    // IG regions: chr 14:106032614-107288051, chr22: 22380474-23265085
    public static void populateImmuneRegions(final List<ChrBaseRegion> regions, final RefGenomeVersion version)
    {
        if(version == RefGenomeVersion.V38)
        {
            regions.add(new ChrBaseRegion("chr2", 88857161, 90315836));
            regions.add(new ChrBaseRegion("chr6", 29722775, 33143325));
            regions.add(new ChrBaseRegion("chr14", 105586437, 106879844));
            regions.add(new ChrBaseRegion("chr22", 22026076, 22922913));
        }
        else
        {
            regions.add(new ChrBaseRegion("2", 89156674, 90538397));
            regions.add(new ChrBaseRegion("6", 29690552, 33111102));
            regions.add(new ChrBaseRegion("14", 106032614, 107288051));
            regions.add(new ChrBaseRegion("22", 22380474, 23265085));
        }
    }



}
