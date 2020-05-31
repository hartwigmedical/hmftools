package com.hartwig.hmftools.linx.ext_compare;

import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.lowerChromosome;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.commons.compress.utils.Lists;

public class CfSvChainData
{
    public final int RearrangeId;
    public final int ChainId;
    public final int[] BreakpointIds;
    public final String[] Chromosomes;
    public final int[] Positions;
    public final byte[] Orientations;
    public final int[] DbBreakendIds;
    public final List<Integer>[] AdjacentBreakendIds;

    private SvVarData mSvData;

    public CfSvChainData(final CfBreakendData[] breakends)
    {
        int startIndex = SE_START;
        int endIndex = SE_END;

        if(breakends[SE_START].Chromosome.equals(breakends[SE_END].Chromosome) && breakends[SE_START].Position > breakends[SE_END].Position
        || lowerChromosome(breakends[SE_END].Chromosome, breakends[SE_START].Chromosome))
        {
            startIndex = SE_END;
            endIndex = SE_START;
        }

        ChainId = breakends[SE_START].ChainId;
        BreakpointIds = new int[] { breakends[startIndex].BreakpointId, breakends[endIndex].BreakpointId};
        RearrangeId = breakends[SE_START].RearrangeId;
        Chromosomes = new String[] { breakends[startIndex].Chromosome, breakends[endIndex].Chromosome};
        Positions = new int[] { breakends[startIndex].Position, breakends[endIndex].Position};
        Orientations = new byte[] { breakends[startIndex].Orientation, breakends[endIndex].Orientation};
        DbBreakendIds = new int[] { breakends[startIndex].DbBreakendId, breakends[endIndex].DbBreakendId };
        AdjacentBreakendIds = new List[] { breakends[startIndex].AdjacentBreakendIds, breakends[endIndex].AdjacentBreakendIds };
        mSvData = null;
    }

    public static final String CF_DATA_DELIMITER = "\t\\s*";

    public final SvVarData getSvData() { return mSvData; }
    public void setSvData(final SvVarData var) { mSvData = var; }

    public boolean matches(final SvVarData var)
    {
        if(var.isSglBreakend())
            return false;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(!var.chromosome(isStart(se)).equals(Chromosomes[se]) || var.position(isStart(se)) != Positions[se]
            || var.orientation(isStart(se)) != Orientations[se])
            {
                return false;
            }
        }

        return true;
    }

    public String toString()
    {
        return String.format("id(%d) chain(%d) pos(%s:%d:%d -> %s:%d:%d) bnds(%d % %d) svId(%d)",
                RearrangeId, ChainId, Chromosomes[SE_START], Orientations[SE_START], Positions[SE_START],
                Chromosomes[SE_END], Orientations[SE_END], Positions[SE_END],
                BreakpointIds[SE_START], BreakpointIds[SE_END], mSvData != null ? mSvData.id() : -1);

    }

}
