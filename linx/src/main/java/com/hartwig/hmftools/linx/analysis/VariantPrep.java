package com.hartwig.hmftools.linx.analysis;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.linx.cn.JcnCalcData;
import com.hartwig.hmftools.linx.cn.SvCNData;
import com.hartwig.hmftools.linx.types.SvVarData;

public class VariantPrep
{
    public static void setSvCopyNumberData(List<SvVarData> svList, final Map<Integer, JcnCalcData> svJcnCalcDataMap,
            final Map<Integer, SvCNData[]> svIdCnDataMap, final Map<String,List<SvCNData>> chrCnDataMap)
    {
        if((svJcnCalcDataMap == null || svJcnCalcDataMap.isEmpty()) && svIdCnDataMap.isEmpty())
            return;

        List<SvCNData> cnDataList = null;
        String currentChromosome = "";
        for(final SvVarData var : svList)
        {
            if(svJcnCalcDataMap != null)
            {
                final JcnCalcData jcnData = svJcnCalcDataMap.get(var.id());
                if (jcnData != null)
                {
                    double estJcn = jcnData.JcnEstimate;
                    double estUncertainty = jcnData.JcnUncertainty;
                    var.setJcnRecalcData(estJcn - estUncertainty, estJcn + estUncertainty);
                }
            }

            final SvCNData[] cnDataPair = svIdCnDataMap.get(var.id());

            if(cnDataPair == null)
                continue;

            for(int be = SE_START; be <= SE_END; ++be)
            {
                if(var.isSglBreakend() && be == SE_END)
                    continue;

                boolean isStart = isStart(be);

                if(!currentChromosome.equals(var.chromosome(isStart)))
                {
                    currentChromosome = var.chromosome(isStart);
                    cnDataList = chrCnDataMap.get(currentChromosome);
                }

                SvCNData cnDataPost = cnDataPair[be];

                if(cnDataList == null || cnDataPost == null)
                    continue;

                SvCNData cnDataPrev = cnDataList.get(cnDataPost.getIndex() - 1);

                var.setCopyNumberData(isStart, cnDataPrev, cnDataPost);
            }
        }
    }



}
