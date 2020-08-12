package com.hartwig.hmftools.linx.ext_compare;

import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.linx.LinxOutput.SUBSET_DELIM;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.sv.SvRegion;

public class AmpliconData
{
    public final List<SvRegion> Regions;
    public final String Type;
    public final int ClusterId;
    public final double MaxCN;
    public final double MaxJCN;
    public final double FoldbackCN;
    public final int Breakends;
    public final int HighBreakends;

    public static final String DATA_DELIM = "\t";

    public static final double MIN_CLUSTER_JCN_THRESHOLD = 3;
    public static final double HIGH_SV_JCN_THRESHOLD = 7;

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
        for(SvRegion region : Regions)
        {
            if(!chromosomes.contains(region.Chromosome))
                chromosomes.add(region.Chromosome);
        }

        return appendStrList(chromosomes, SUBSET_DELIM);
    }

    public static final String extractSampleId(final Map<String,Integer> fieldsIndexMap, final String[] items)
    {
        String sampleId = items[fieldsIndexMap.get("pair")];

        if(sampleId.endsWith(".2"))
            return sampleId.substring(0, sampleId.length() - 2) + "TII";
        else if(sampleId.endsWith(".3"))
            return sampleId.substring(0, sampleId.length() - 2) + "TIII";
        else if(sampleId.endsWith(".4"))
            return sampleId.substring(0, sampleId.length() - 2) + "TIV";
        else
            return sampleId + "T";
    }

    public static AmpliconData from(final Map<String,Integer> fieldsIndexMap, final String[] items)
    {
        //pair	seqnames	start	end	strand	width	type	footprint	ev.id	id	max.cn	max.jcn	cluster
        final String chromosome = items[fieldsIndexMap.get("seqnames")];
        final int[] positions = { Integer.parseInt(items[fieldsIndexMap.get("start")]), Integer.parseInt(items[fieldsIndexMap.get("end")]) };
        final SvRegion region = new SvRegion(chromosome, positions);
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
