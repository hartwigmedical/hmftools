package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.StructuralVariant;

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

    public void populateAssemblyLinks(final GermlineSV germlineSV)
    {
        final StructuralVariant sv = germlineSV.sv();
        final AssemblyData asmData = mSvAssemblyData.stream().filter(x -> x.VcfId.equals(sv.id())).findFirst().orElse(null);

        if(asmData == null)
            return;

        germlineSV.setAssemblySvId(SE_START, asmData.getLinkedSvIds()[SE_START]);
        germlineSV.setAssemblySvId(SE_END, asmData.getLinkedSvIds()[SE_END]);
    }

}
