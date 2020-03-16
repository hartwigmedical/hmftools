package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.seIndex;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

public class AssemblyData
{
    public final String VcfId;
    public final List<List<String>> BEID;
    public final List<List<String>> BEIDL;

    private String[] mLinkedSvId;

    public AssemblyData(final String vcfId, final String[] beIdStr, final String[] beIdlStr)
    {
        VcfId = vcfId;

        BEID = Lists.newArrayListWithCapacity(SE_PAIR);
        BEIDL = Lists.newArrayListWithCapacity(SE_PAIR);

        for(int se = SE_START; se <= SE_END; ++se)
        {
            final String[] beIdData = beIdStr[se]
                    .replaceAll("\\[", "")
                    .replaceAll("]", "")
                    .split(",");

            List<String> beIds = Lists.newArrayList();

            if(beIdData.length > 1 || !beIdData[0].equals(""))
                Arrays.stream(beIdData).forEach(x -> beIds.add(x));

            BEID.add(beIds);

            final String[] beIdlData = beIdlStr[se]
                    .replaceAll("\\[", "")
                    .replaceAll("]", "")
                    .split(",");

            List<String> beIdls = Lists.newArrayList();

            if(beIdlData.length > 1 || !beIdlData[0].equals(""))
                Arrays.stream(beIdlData).forEach(x -> beIdls.add(x));

            BEIDL.add(beIdls);
        }

        mLinkedSvId = new String[SE_PAIR];
        mLinkedSvId[SE_START] = "";
        mLinkedSvId[SE_END] = "";
    }

    public void setLinkedData(final String svId, boolean linkedOnStart)
    {
        mLinkedSvId[seIndex(linkedOnStart)] = svId;
    }

    public final String[] getLinkedSvIds() { return mLinkedSvId; }

    public boolean hasMatch(final AssemblyData other, int seIndex, int otherSeIndex)
    {
        if(other.VcfId.equals(VcfId))
            return false;

        boolean hasBeIdMatch = false;
        for(final String beId : BEID.get(seIndex))
        {
            if(other.BEID.get(otherSeIndex).stream().anyMatch(x -> x.equals(beId)))
            {
                hasBeIdMatch = true;
                break;
            }
        }

        if(!hasBeIdMatch)
            return false;

        boolean hasBeIdlMatch = false;
        for(final String beIdl : BEIDL.get(seIndex))
        {
            if(other.BEIDL.get(otherSeIndex).stream().anyMatch(x -> x.equals(beIdl)))
            {
                hasBeIdlMatch = true;
                break;
            }
        }

        return hasBeIdlMatch;
    }

    public static void annotateAssembledLinks(final List<AssemblyData> svAssemblyData)
    {
        for(int i = 0; i < svAssemblyData.size() - 1; ++i)
        {
            AssemblyData asmData1 = svAssemblyData.get(i);
            boolean linked = false;

            for(int j = i+1; j < svAssemblyData.size(); ++j)
            {
                AssemblyData asmData2 = svAssemblyData.get(j);

                for(int se1 = SE_START; se1 <= SE_END; ++se1)
                {
                    for(int se2 = SE_START; se2 <= SE_END; ++se2)
                    {
                        if(asmData1.hasMatch(asmData2, se1, se2))
                        {
                            asmData1.setLinkedData(asmData2.VcfId, isStart(se1));
                            asmData2.setLinkedData(asmData1.VcfId, isStart(se2));
                            linked = true;
                            break;
                        }
                    }

                    if(linked)
                        break;
                }

                if(linked)
                    break;
            }
        }
    }

    public static void populateAssemblyLinks(final List<AssemblyData> svAssemblyData, final GermlineSV germlineSV)
    {
        final StructuralVariant sv = germlineSV.sv();
        final AssemblyData asmData = svAssemblyData.stream().filter(x -> x.VcfId.equals(sv.id())).findFirst().orElse(null);

        if(asmData == null)
            return;

        germlineSV.setAssemblySvId(SE_START, asmData.getLinkedSvIds()[SE_START]);
        germlineSV.setAssemblySvId(SE_END, asmData.getLinkedSvIds()[SE_END]);
    }
}
