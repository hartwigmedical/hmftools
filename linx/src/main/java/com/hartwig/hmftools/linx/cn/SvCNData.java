package com.hartwig.hmftools.linx.cn;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.linx.types.SvVarData;

public class SvCNData {

    final private int mId;
    final public String Chromosome;
    final public int StartPos;
    final public int EndPos;
    final public double CopyNumber;
    final public String SegStart;
    final public String SegEnd;
    final public String Method;
    final public int BafCount;
    final public double ObservedBaf;
    final public double ActualBaf;
    final public int DepthWindowCount;

    private int mIndex; // in the copy number table

    private StructuralVariantData mSvData; // linked if known
    private boolean mSvLinkOnStart;

    public SvCNData(int id, final String chromosome, final int startPos, final int endPos,
            final double copyNumber, final String segStart,  final String segEnd,  final int bafCount,
            final double actualBaf,  final int depthWindowCount)
    {
        mId = id;
        Chromosome = chromosome;
        StartPos = startPos;
        EndPos = endPos;
        CopyNumber = copyNumber;
        SegStart = segStart;
        SegEnd = segEnd;
        Method = CopyNumberMethod.STRUCTURAL_VARIANT.toString();
        BafCount = bafCount;
        ObservedBaf = 0;
        ActualBaf = actualBaf;
        DepthWindowCount = depthWindowCount;
        mIndex = 0;
    }

    public SvCNData(final PurpleCopyNumber record, int id)
    {
        mId = id;
        Chromosome = record.chromosome();
        StartPos = (int)record.start();
        EndPos = (int)record.end();
        CopyNumber = record.averageTumorCopyNumber();
        SegStart = record.segmentStartSupport().toString();
        SegEnd = record.segmentEndSupport().toString();
        Method = record.method().toString();
        BafCount = record.bafCount();
        ObservedBaf = record.averageObservedBAF();
        ActualBaf = record.averageActualBAF();
        DepthWindowCount = record.depthWindowCount();
        mIndex = 0;
    }

    public int id() { return mId; }
    public int position(boolean useStart) { return useStart ? StartPos : EndPos; }
    public int length() { return EndPos - StartPos; }

    public int getIndex() { return mIndex; }
    public void setIndex(int index) { mIndex = index; }

    public double majorAlleleJcn() { return ActualBaf * CopyNumber; }
    public double minorAlleleJcn() { return max((1 - ActualBaf) * CopyNumber,0); }

    public void setStructuralVariantData(final StructuralVariantData svData, boolean linkOnStart)
    {
        mSvData = svData;
        mSvLinkOnStart = linkOnStart;
    }

    public final StructuralVariantData getStructuralVariantData() { return mSvData; }
    public boolean svLinkOnStart() { return mSvLinkOnStart; }

    public boolean matchesSegment(SegmentSupport segment, boolean isStart)
    {
        return isStart ? SegStart.equals(segment.toString()) : SegEnd.equals(segment.toString());
    }

    public static boolean isSvSegment(final SegmentSupport segment)
    {
        return (segment == SegmentSupport.BND || segment == SegmentSupport.INV || segment == SegmentSupport.INS
                || segment == SegmentSupport.DEL || segment == SegmentSupport.DUP || segment == SegmentSupport.SGL
                || segment == SegmentSupport.MULTIPLE);
     }

    public boolean matchesSV(boolean isStart)
    {
        return isStart ? isSvSegment(SegmentSupport.valueOf(SegStart)) : isSvSegment(SegmentSupport.valueOf(SegEnd));
    }

    public final String toString()
    {
        return String.format("chr(%s) pos(%d - %d) segs(%s - %s) cn(%.2f) baf(%.2f)",
                Chromosome, StartPos, EndPos, SegStart, SegEnd, CopyNumber, ActualBaf);
    }

    public static void setSvCopyNumberData(
            final List<SvVarData> svList, final Map<Integer,JcnCalcData> svJcnCalcDataMap,
            final Map<Integer,SvCNData[]> svIdCnDataMap, final Map<String,List<SvCNData>> chrCnDataMap)
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
