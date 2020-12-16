package com.hartwig.hmftools.linx.ext_compare;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.lowerChromosome;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.linx.ext_compare.CfBreakendData.NO_ID;
import static com.hartwig.hmftools.linx.ext_compare.CfDbMatchType.CF_ONLY;
import static com.hartwig.hmftools.linx.ext_compare.CfDbMatchType.DIFF;
import static com.hartwig.hmftools.linx.ext_compare.CfDbMatchType.LINX_ONLY;
import static com.hartwig.hmftools.linx.ext_compare.CfDbMatchType.MATCH;
import static com.hartwig.hmftools.linx.types.LinxConstants.NO_DB_MARKER;

import java.util.List;

import com.hartwig.hmftools.common.utils.sv.StartEndPair;
import com.hartwig.hmftools.linx.types.DbPair;
import com.hartwig.hmftools.linx.types.SvVarData;

public class CfSvData
{
    public final int RearrangeId;
    public final int ChainId;
    public final int[] BreakpointIds;
    public final String[] Chromosomes;
    public final int[] Positions;
    public final byte[] Orientations;
    public final int[] DbBreakendIds;
    public final StartEndPair<List<Integer>> AdjacentBreakendIds;

    private SvVarData mSvData;
    private final SvVarData[] mDbSvData;
    private final int[] mDbLengths;

    private final CfLink[] mChainLinks;

    public CfSvData(final CfBreakendData[] breakends)
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
        AdjacentBreakendIds = new StartEndPair<>(breakends[startIndex].AdjacentBreakendIds, breakends[endIndex].AdjacentBreakendIds);

        mSvData = null;
        mDbSvData = new SvVarData[] { null, null };
        mDbLengths = new int[] { NO_DB_MARKER, NO_DB_MARKER };
        mChainLinks = new CfLink[] { null, null };
    }

    public static final String CF_DATA_DELIMITER = "\t\\s*";

    public final SvVarData getSvData() { return mSvData; }
    public void setSvData(final SvVarData var) { mSvData = var; }

    public final SvVarData[] getDbSvData() { return mDbSvData; }
    public final int[] getDbLengths() { return mDbLengths; }

    public final CfLink[] getChainLinks() { return mChainLinks; }
    public void setChainLink(final CfLink link, int index) { mChainLinks[index] = link; }

    public void setDbSvData(int se, final SvVarData var, int length)
    {
        mDbSvData[se] = var;
        mDbLengths[se] = length;
    }

    public boolean svMatched() { return mSvData != null; }

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

    public int calcDbLength(final CfSvData other, int se, int seOther)
    {
        if(Orientations[se] == other.Orientations[seOther])
            return NO_DB_MARKER;

        if(Orientations[se] == POS_ORIENT)
            return other.Positions[seOther] - Positions[se];
        else
            return Positions[se] - other.Positions[seOther];
    }

    public CfDbMatchType[] getDeletionBridgeMatchTypes()
    {
        final CfDbMatchType[] dbMatchTypes = new CfDbMatchType[] { CfDbMatchType.NONE, CfDbMatchType.NONE };

        for(int se = SE_START; se <= SE_END; ++se)
        {
            final DbPair dbLink = mSvData.getDBLink(isStart(se));

            if(dbLink == null)
            {
                if(DbBreakendIds[se] == NO_ID)
                    continue;

                dbMatchTypes[se] = CF_ONLY;
            }
            else if(DbBreakendIds[se] == NO_ID)
            {
                dbMatchTypes[se] = LINX_ONLY;
            }
            else
            {
                if(dbLink.getOtherSV(mSvData) == mDbSvData[se])
                {
                    dbMatchTypes[se] = MATCH;
                }
                else
                {
                    dbMatchTypes[se] = DIFF;
                }
            }
        }

        return dbMatchTypes;
    }

    public String toString()
    {
        return String.format("id(%d) chain(%d) pos(%s:%d:%d -> %s:%d:%d) bnds(%d & %d) svId(%d)",
                RearrangeId, ChainId, Chromosomes[SE_START], Orientations[SE_START], Positions[SE_START],
                Chromosomes[SE_END], Orientations[SE_END], Positions[SE_END],
                BreakpointIds[SE_START], BreakpointIds[SE_END], mSvData != null ? mSvData.id() : -1);

    }

}
