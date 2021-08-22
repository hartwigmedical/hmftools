package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.AS;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.BEID;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.BEIDL;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.CAS;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.RAS;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.StructuralVariant;

import htsjdk.variant.variantcontext.CommonInfo;

public class LinkAnalyser
{
    private final List<AssemblyData> mSvAssemblyData;

    public LinkAnalyser()
    {
        mSvAssemblyData = Lists.newArrayList();

    }

    public void clear()
    {
        mSvAssemblyData.clear();
    }

    public void add(final AssemblyData asmData)
    {
        mSvAssemblyData.add(asmData);
    }

    public void cacheAssemblyData(final SvData sv)
    {
        final CommonInfo variantCI = sv.contextStart().getCommonInfo();

        if(variantCI.getAttributeAsInt(AS, 0) + variantCI.getAttributeAsInt(CAS, 0)
                + variantCI.getAttributeAsInt(RAS, 0) < 2)
        {
            return;
        }

        // cache assembly info
        final String startBeId = sv.contextStart().getAttributeAsString(BEID,"");
        final String startBeIdl = sv.contextStart().getAttributeAsString(BEIDL,"");

        final String endBeId = sv.contextEnd() != null ? sv.contextEnd().getAttributeAsString(BEID,"") : "";
        final String endBeIdl = sv.contextEnd() != null ? sv.contextEnd().getAttributeAsString(BEIDL,"") : "";

        if(startBeId.isEmpty() && endBeId.isEmpty())
            return;

        final String[] beIdStr = {startBeId, endBeId};
        final String[] beIdlStr = {startBeIdl, endBeIdl};

        AssemblyData asmData = new AssemblyData(sv.id(), beIdStr, beIdlStr);
        add(asmData);
    }


    public void annotateAssembledLinks()
    {
        for(int i = 0; i < mSvAssemblyData.size() - 1; ++i)
        {
            AssemblyData asmData1 = mSvAssemblyData.get(i);
            boolean linked = false;

            for(int j = i+1; j < mSvAssemblyData.size(); ++j)
            {
                AssemblyData asmData2 = mSvAssemblyData.get(j);

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

    public void populateAssemblyLinks(final SvData svData)
    {
        final AssemblyData asmData = mSvAssemblyData.stream().filter(x -> x.VcfId.equals(svData.id())).findFirst().orElse(null);

        if(asmData == null)
            return;

        svData.setAssemblySvId(SE_START, asmData.getLinkedSvIds()[SE_START]);
        svData.setAssemblySvId(SE_END, asmData.getLinkedSvIds()[SE_END]);
    }

}
