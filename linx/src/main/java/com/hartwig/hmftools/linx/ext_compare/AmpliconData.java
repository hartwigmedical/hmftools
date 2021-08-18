package com.hartwig.hmftools.linx.ext_compare;

import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.linx.LinxOutput.ITEM_DELIM_CHR;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class AmpliconData
{
    public final List<ChrBaseRegion> Regions;
    public final String Type;
    public final int ClusterId;
    public final double MaxCN;
    public final double MaxJCN;
    public final double FoldbackCN;
    public final int Breakends;
    public final int HighBreakends;

    public static final String JABBA_DATA_DELIM = "\t";
    public static final String AA_DATA_DELIM = ",";

    public static final String AMPLICON_SOURCE_AMP_ARCHITECT = "AMPLICON_ARCHITECT";
    public static final String AMPLICON_SOURCE_JABBA = "JABBA";

    public AmpliconData(final int clusterId, final String type, final double maxCN, final double maxJCN,
            final double foldbackCN, final int breakends, final int highBreakends)
    {
        Regions = Lists.newArrayList();
        Type = type;
        ClusterId = clusterId;
        MaxCN = maxCN;
        MaxJCN = maxJCN;
        FoldbackCN = foldbackCN;
        Breakends = breakends;
        HighBreakends = highBreakends;
    }

    public String chromosomes()
    {
        final List<String> chromosomes = Lists.newArrayList();
        for(ChrBaseRegion region : Regions)
        {
            if(!chromosomes.contains(region.Chromosome))
                chromosomes.add(region.Chromosome);
        }

        return appendStrList(chromosomes, ITEM_DELIM_CHR);
    }

    public static String extractSampleId(final String dataSource, final Map<String,Integer> fieldsIndexMap, final String[] items)
    {
        if(dataSource.equals(AMPLICON_SOURCE_AMP_ARCHITECT))
        {
            String sampleId = items[fieldsIndexMap.get("SampleId")];

            return sampleId.endsWith("T") ? sampleId : sampleId + "T";
        }

        String sampleId = items[fieldsIndexMap.get("pair")];

        if(sampleId.endsWith("T"))
            return sampleId;

        if(sampleId.endsWith(".2"))
            return sampleId.substring(0, sampleId.length() - 2) + "TII";
        else if(sampleId.endsWith(".3"))
            return sampleId.substring(0, sampleId.length() - 2) + "TIII";
        else if(sampleId.endsWith(".4"))
            return sampleId.substring(0, sampleId.length() - 2) + "TIV";
        else
            return sampleId + "T";
    }

    public static AmpliconData from(final String dataSource, final Map<String,Integer> fieldsIndexMap, final String[] items)
    {
        if(dataSource.equals(AMPLICON_SOURCE_AMP_ARCHITECT))
        {
            if(!fieldsIndexMap.containsKey("SampleId") || !fieldsIndexMap.containsKey("AmpliconIntervals")
            || !fieldsIndexMap.containsKey("AmpliconClusterId") || !fieldsIndexMap.containsKey("AmpliconClassification"))
            {
                return null;
            }

            // SampleBarcode,SampleId,AmpliconClusterId,AmpliconClassification,AmpliconIntervals

            final String type = items[fieldsIndexMap.get("AmpliconClassification")];
            final int clusterId = Integer.parseInt(items[fieldsIndexMap.get("AmpliconClusterId")]);

            AmpliconData ampData = new AmpliconData(clusterId, type, 0, 0, 0, 0, 0);

            // intervals: 3:174885846-197847500;11:40144508-40344836
            final String[] ampRegions = items[fieldsIndexMap.get("AmpliconIntervals")].split(";", -1);

            for(String ampRegion : ampRegions)
            {
                final String[] regionItems = ampRegion.split(":");

                if(regionItems.length != 3)
                    return null;

                final String chromosome = regionItems[0];
                final int[] positions = { Integer.parseInt(regionItems[1]), Integer.parseInt(regionItems[2]) };
                ampData.Regions.add(new ChrBaseRegion(chromosome, positions));
            }

            return ampData;
        }

        //pair	seqnames	start	end	strand	width	type	footprint	ev.id	id	max.cn	max.jcn	cluster
        final String chromosome = items[fieldsIndexMap.get("seqnames")];
        final int[] positions = { Integer.parseInt(items[fieldsIndexMap.get("start")]), Integer.parseInt(items[fieldsIndexMap.get("end")]) };
        final ChrBaseRegion region = new ChrBaseRegion(chromosome, positions);
        final String type = items[fieldsIndexMap.get("type")];
        final int clusterId = Integer.parseInt(items[fieldsIndexMap.get("cluster")]);
        final double maxCN = Integer.parseInt(items[fieldsIndexMap.get("max.cn")]);
        final double maxJCN = Integer.parseInt(items[fieldsIndexMap.get("max.jcn")]);
        final double foldbackCN = Integer.parseInt(items[fieldsIndexMap.get("fbi.cn")]);
        final int breakends = Integer.parseInt(items[fieldsIndexMap.get("n.jun")]);
        final int highBreakends = Integer.parseInt(items[fieldsIndexMap.get("n.jun.high")]);

        AmpliconData ampData = new AmpliconData(clusterId, type, maxCN, maxJCN, foldbackCN, breakends, highBreakends);
        ampData.Regions.add(region);
        return ampData;
    }

    public String toString()
    {
        return String.format("%d: %s regions(%d) cn(%.1f) jcn(%.1f) breakends(%d high=%d)",
                ClusterId, Type, Regions.size(), MaxCN, MaxJCN, Breakends, HighBreakends);
    }



}
